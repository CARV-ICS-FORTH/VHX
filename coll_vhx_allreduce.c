

#include "ompi_config.h"
#include "mpi.h"

#include "ompi/constants.h"
#include "ompi/datatype/ompi_datatype.h"
#include "ompi/communicator/communicator.h"
#include "ompi/op/op.h"

#include "ompi/mca/coll/base/base.h"

#include "opal/mca/rcache/base/base.h"
#include "opal/util/minmax.h"

#include "coll_vhx.h"



// -----------------------------

int mca_coll_vhx_allreduce_internal(const void * sbuf, void * rbuf,
  int count, ompi_datatype_t * datatype, ompi_op_t * op,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module, int bcast) {

  int rank = ompi_comm_rank(ompi_comm);
  int comm_size = ompi_comm_size(ompi_comm);
  size_t dtype_size, bytes_total;
  ompi_datatype_type_size(datatype, & dtype_size);
  bytes_total = count * dtype_size;
  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  if (((mca_coll_vhx_module_t * ) module) -> initialized == false) {
    int ret = vhx_init(module, ompi_comm);
    if (ret != 0) return OMPI_ERROR;
  }
  int pvt_seq;;;
  int hier_size = vhx_module -> hierarchy_size;
  int first_reduction_done = 0;
  int do_cico = (bytes_total <= OMPI_vhx_CICO_MAX);

  for (int i = 0; i < hier_size; i++) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
    hier_group -> shared_ctrl_vars[rank].rbuf_vaddr = rbuf;
    hier_group -> shared_ctrl_vars[rank].sbuf_vaddr = sbuf;

    if (rank == hier_group -> leader) {
      pvt_seq = ++(vhx_module -> hier_groups[i].shared_ctrl_vars[rank].coll_seq);

      for (int j = 0; j < hier_group -> real_size; j++) { //waitng for SEQ wave

        if (hier_group -> real_members[j] == rank)
          continue;

        while (hier_group -> shared_ctrl_vars[hier_group -> real_members[j]].coll_seq != pvt_seq);
        opal_atomic_rmb();
        if (!do_cico) {
          mca_coll_vhx_get_rank_reg(hier_group -> real_members[j], hier_group -> shared_ctrl_vars[hier_group -> real_members[j]].sbuf_vaddr,
            bytes_total, & (hier_group -> sbuf_regs[hier_group -> real_members[j]]), vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[hier_group -> real_members[j]]);
          mca_coll_vhx_get_rank_reg(hier_group -> real_members[j], hier_group -> shared_ctrl_vars[hier_group -> real_members[j]].rbuf_vaddr,
            bytes_total, & (hier_group -> rbuf_regs[hier_group -> real_members[j]]), vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[hier_group -> real_members[j]]);
        }

        if (!first_reduction_done) {
          ompi_3buff_op_reduce(op, do_cico ? vhx_module -> neighbour_cico_buffers[hier_group -> real_members[j]] : hier_group -> neighbour_sbufs[hier_group -> real_members[j]], sbuf, (do_cico) ? vhx_module -> cico_buffer : rbuf, count, datatype);
          first_reduction_done++;
        } else if (i == 0) //in the first level we use sendbufs as input in no CICO
          ompi_op_reduce(op, do_cico ? vhx_module -> neighbour_cico_buffers[hier_group -> real_members[j]] : hier_group -> neighbour_sbufs[hier_group -> real_members[j]], (do_cico) ? vhx_module -> cico_buffer : rbuf, count, datatype);
        else // in the next levels, rnuf as eused as input sincethey contain the product of previous reductions in no CICO
          ompi_op_reduce(op, do_cico ? vhx_module -> neighbour_cico_buffers[hier_group -> real_members[j]] : hier_group -> neighbour_rbufs[hier_group -> real_members[j]], (do_cico) ? vhx_module -> cico_buffer : rbuf, count, datatype);

      }
      if (i == hier_size - 1 && do_cico)
        memcpy(rbuf, vhx_module -> cico_buffer, bytes_total); // in no CICO we need to copy the product of the final reduction to its recv buffer
      opal_atomic_wmb();
      for (int j = 0; j < hier_group -> real_size; j++)
        hier_group -> shared_ctrl_vars[hier_group -> real_members[j]].coll_ack = pvt_seq;

    } else if (hier_group -> real_members_bitmap[rank] == 1) {
      if (do_cico && i == 0)
        memcpy(vhx_module -> cico_buffer, sbuf, bytes_total);
      opal_atomic_wmb();

      pvt_seq = ++(vhx_module -> hier_groups[i].shared_ctrl_vars[rank].coll_seq);

      while (hier_group -> shared_ctrl_vars[rank].coll_ack != pvt_seq);

    }
  }

  if (bcast) {

    pvt_seq = (vhx_module -> hier_groups[0].shared_ctrl_vars[rank].coll_seq) + 1;

    int root = vhx_module -> hier_groups[hier_size - 1].leader;
    for (int i = hier_size - 1; i >= 0; i--) {
      vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

      (hier_group -> shared_ctrl_vars[rank].coll_seq) = pvt_seq;
      if (rank == hier_group -> leader) {

        opal_atomic_wmb();

        vhx_module -> hier_groups[i].shared_ctrl_vars[rank].coll_ack = pvt_seq; // we need the root's ack counter to be equal to pvt for consistency
        if (i == 0 && rank != root && do_cico)
          memcpy(rbuf, (char * )(vhx_module -> cico_buffer), bytes_total);

      } else if (hier_group -> real_members_bitmap[rank] == 1) {

        WAIT_FLAG( & (vhx_module -> hier_groups[i].shared_ctrl_vars[hier_group -> leader].coll_seq), pvt_seq);
        opal_atomic_rmb();

        if (do_cico)
          memcpy((i == 0) ? rbuf : (char * )(vhx_module -> cico_buffer), (char * )(vhx_module -> neighbour_cico_buffers[hier_group -> leader]), bytes_total); // on the bottom we need to write to buf even in cico scenarios
        else {

          mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> shared_ctrl_vars[hier_group -> leader].rbuf_vaddr,
            bytes_total, & (hier_group -> rbuf_regs[hier_group -> leader]), vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[hier_group -> leader]);

          memcpy(rbuf, (char * )(hier_group -> neighbour_rbufs[hier_group -> leader]), bytes_total);

        }
        opal_atomic_wmb();

      }
    }

    for (int i = 0; i < hier_size; i++) {
      vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
      if (rank == hier_group -> leader)
        for (int j = 0; j < hier_group -> real_size; j++) { //Waiting for the acks of the hierarchy group
          if (hier_group -> real_members[j] == rank)
            continue;
          WAIT_FLAG( & (hier_group -> shared_ctrl_vars[hier_group -> real_members[j]].coll_ack), pvt_seq);
        }
      hier_group -> shared_ctrl_vars[rank].coll_ack = pvt_seq; // we need the root's ack counter to be equal to pvt for consistency
    }
  }

  return OMPI_SUCCESS;
}

int mca_coll_vhx_allreduce(const void * sbuf, void * rbuf,
  int count, ompi_datatype_t * datatype, ompi_op_t * op,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module) {

  return mca_coll_vhx_allreduce_internal(sbuf, rbuf,
    count, datatype, op, ompi_comm, module, true);
}
