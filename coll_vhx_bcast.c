

#include "ompi_config.h"
#include "mpi.h"

#include "ompi/constants.h"
#include "ompi/datatype/ompi_datatype.h"
#include "ompi/communicator/communicator.h"
#include "opal/util/show_help.h"
#include "opal/util/minmax.h"

#include "coll_vhx.h"

//the root becomes the leader in all hierarchy groups it participates
int set_leader(int root, mca_coll_base_module_t * module, ompi_communicator_t * ompi_comm) {

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i < hier_size; i++) {
	vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

    if (hier_group -> members_bitmap[root] == 1 && hier_group -> real_members_bitmap[root] == 0) {
      hier_group -> real_members_bitmap[root] = 1;

      hier_group -> leader = root;
    }
	
	
  }
  return OMPI_SUCCESS;
}


		
    
	
	
  

int set_vaddr(int my_rank, mca_coll_base_module_t * module, void * sbuf){
	
	vhx_module_t * vhx_module = (vhx_module_t * ) module;
	int hier_size = vhx_module -> hierarchy_size;
	
	for (int i = 0; i < hier_size; i++) {
		vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

		
		hier_group->shared_ctrl_vars[my_rank].sbuf_vaddr = sbuf;
    }
	
	
}

int mca_coll_vhx_bcast(void * buf, int count, ompi_datatype_t * datatype, int root,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module) {

  int rank = ompi_comm_rank(ompi_comm);
  int comm_size = ompi_comm_size(ompi_comm);

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  if (((mca_coll_vhx_module_t * ) module) -> initialized == false) {
    int ret = vhx_init(module, ompi_comm);
    if (ret != 0) return OMPI_ERROR;
  }
  int hier_size = vhx_module -> hierarchy_size;

  int pvt_seq = (vhx_module -> hier_groups[0].shared_ctrl_vars[rank].coll_seq) + 1;
  size_t dtype_size, bytes_total;
  ompi_datatype_type_size(datatype, & dtype_size);
  bytes_total = count * dtype_size;

  bool do_cico = (bytes_total <= OMPI_vhx_CICO_MAX);
  //do_cico = 1;	//TODO, XPMEM support temporarily broken due to recent "refactoring" after discovering erroneous execution, will get fixed soon

  set_leader(root, vhx_module, ompi_comm);
  if (do_cico && root == rank)
    memcpy((char * )(vhx_module -> cico_buffer), buf, bytes_total);
  if (!do_cico)
	set_vaddr(rank, module, buf);

  for (int i = hier_size - 1; i >= 0; i--) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

    
      (hier_group -> shared_ctrl_vars[rank].coll_seq) = pvt_seq;
    if (rank == hier_group -> leader) {

      opal_atomic_wmb();

      vhx_module -> hier_groups[i].shared_ctrl_vars[rank].coll_ack = pvt_seq; // we need the root's ack counter to be equal to pvt for consistency
      if (i == 0 && rank != root && do_cico)
        memcpy(buf, (char * )(vhx_module -> cico_buffer), bytes_total);

      for (int i = 0; i < hier_group -> real_size; i++)  //Waiting for the acks of the hierarchy group
        WAIT_FLAG( & (hier_group -> shared_ctrl_vars[hier_group -> real_members[i]].coll_ack), pvt_seq);

    } else if (hier_group -> real_members_bitmap[rank] == 1) {

      WAIT_FLAG( & (vhx_module -> hier_groups[i].shared_ctrl_vars[hier_group -> leader].coll_seq), pvt_seq);
	  opal_atomic_rmb();
 
	  if (do_cico)
		memcpy((i == 0) ? buf : (char * )(vhx_module -> cico_buffer), (char * )(vhx_module -> neighbour_cico_buffers[hier_group -> leader]), bytes_total); // on the bottom we need to write to buf even in cico scenarios
      else{
		
		mca_coll_vhx_get_rank_reg(hier_group->leader, hier_group->shared_ctrl_vars[hier_group->leader].sbuf_vaddr,
				bytes_total, &(hier_group->sbuf_regs[hier_group->leader]), vhx_module, ompi_comm, &hier_group->neighbour_sbufs[hier_group->leader]);  
	
		 memcpy( buf, (char * )(hier_group->neighbour_sbufs[hier_group->leader]), bytes_total); 
		 
	  }	opal_atomic_wmb();
      hier_group -> shared_ctrl_vars[rank].coll_ack = pvt_seq;

    }
  }

  return OMPI_SUCCESS;

}
