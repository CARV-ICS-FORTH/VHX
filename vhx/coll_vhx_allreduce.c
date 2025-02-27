
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


void print_memory_decimal(const void * mem, size_t n) {
  float * byte = (float * ) mem;

  for (size_t i = 0; i < n / sizeof(float); ++i) {
    printf("%f ", byte[i]);
  }

  printf("\n");
}


int init_bytes_ready_allreduce(int my_rank, vhx_module_t * module, size_t chunks, int pvt_seq) {

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i < hier_size; i++) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
		  WAIT_FLAG(& hier_group -> shared_ctrl_vars[0].coll_ack, pvt_seq - 1);

    if (hier_group -> real_members_bitmap[my_rank] == 1) {

      volatile size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = & (hier_group -> shared_ctrl_vars[0].bytes_available);
      __atomic_store_n(tmp_bytes_available, 0, __ATOMIC_RELAXED);

      if (i == 0) {
        __atomic_store_n( & (hier_group -> shared_ctrl_vars[my_rank].reduce_ready_chunks), chunks, __ATOMIC_RELAXED);

      } else
        __atomic_store_n( & (hier_group -> shared_ctrl_vars[my_rank].reduce_ready_chunks), 0, __ATOMIC_RELAXED);
      if (hier_group -> leader != my_rank)
        break;
    }
  }
  return 0;

}

int set_bytes_ready_allreduce(int my_rank, vhx_hier_group_t * hier_group, size_t chunks) {

  __atomic_store_n( & (hier_group -> shared_ctrl_vars[my_rank].reduce_ready_chunks), chunks, __ATOMIC_RELAXED);

}

int get_left_neighbour(vhx_hier_group_t * hier_group, int rank) {
  int i = 1;
  while (1)
    if (hier_group -> real_members_bitmap[rank - i] != 1)
      i++;
    else
      break;
  return rank - i;

}
// -----------------------------

int reduce_internal_cico(const void * sbuf, void * rbuf,
  int count, ompi_datatype_t * datatype, ompi_op_t * op,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module, int pvt_seq) {
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
  int hier_size = vhx_module -> hierarchy_size;
  int do_cico = (bytes_total <= OMPI_vhx_CICO_MAX);
  size_t chunks = (OMPI_vhx_CHUNK_SIZE && OMPI_vhx_CHUNK_SIZE < bytes_total) ? bytes_total / OMPI_vhx_CHUNK_SIZE : 1;
  int remainder = (OMPI_vhx_CHUNK_SIZE && OMPI_vhx_CHUNK_SIZE < bytes_total) ? bytes_total % OMPI_vhx_CHUNK_SIZE : 0;
  int first_reduction_done = -chunks + 1;

  int has_remainder = remainder ? 1 : 0;
  init_bytes_ready_allreduce(rank, vhx_module, chunks, pvt_seq);
  for (int i = 0; i < hier_size; i++) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
    hier_group -> shared_ctrl_vars[rank].rbuf_vaddr = rbuf;
    hier_group -> shared_ctrl_vars[rank].sbuf_vaddr = (void*) sbuf;
    if (rank == hier_group -> leader) {
	//	printf(" i am leader %d lvl %d\n", rank, i);
      hier_group -> shared_ctrl_vars[rank].coll_seq = pvt_seq;
      for (int j = 0; j < hier_group -> real_size; j++) { //waitng for SEQ wave

        if (hier_group -> real_members[j] == rank)
          continue;
	//	printf(" i am leader %d lvl %d  pre hhh\n", rank, i);
        while (hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_seq < pvt_seq);
        opal_atomic_rmb();
	//	printf(" i am leader %d lvl %d  hhh\n", rank, i);

        if (!do_cico) {
          mca_coll_vhx_get_rank_reg(hier_group -> real_members[j], hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].sbuf_vaddr,
            bytes_total, & (hier_group -> sbuf_regs[hier_group -> real_members[j]]), (mca_coll_base_module_t *)vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[hier_group -> real_members[j]]);
          mca_coll_vhx_get_rank_reg(hier_group -> real_members[j], hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].rbuf_vaddr,
            bytes_total, & (hier_group -> rbuf_regs[hier_group -> real_members[j]]), (mca_coll_base_module_t *)vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[hier_group -> real_members[j]]);
        }

        int reduced_chunks = 0;
        while (reduced_chunks < chunks + has_remainder) {
          size_t offset = 0;
          int reduce_count = count;
          int chunks_ready = __atomic_load_n( & (hier_group -> shared_ctrl_vars[hier_group -> real_members[j]].reduce_ready_chunks), __ATOMIC_RELAXED);
          opal_atomic_rmb();
          if (chunks_ready <= reduced_chunks)
            continue;
          if (OMPI_vhx_CHUNK_SIZE) {
            offset = reduced_chunks * (OMPI_vhx_CHUNK_SIZE);

            reduce_count = (reduced_chunks == chunks && has_remainder) ? remainder : opal_min(OMPI_vhx_CHUNK_SIZE, bytes_total) / dtype_size;


          }
          if (first_reduction_done != 1) {
            ompi_3buff_op_reduce(op, do_cico ? (char * ) vhx_module -> neighbour_cico_buffers[hier_group -> real_members[j]] + offset : (char * ) hier_group -> neighbour_sbufs[hier_group -> real_members[j]] + offset, (char * ) sbuf + offset, (do_cico) ? (char * ) vhx_module -> cico_buffer + offset : (char * ) rbuf + offset, reduce_count, datatype);
            first_reduction_done++;
          } else if (i == 0) //in the first level we use sendbufs as input in no CICO
            ompi_op_reduce(op, do_cico ? (char * )(vhx_module -> neighbour_cico_buffers[hier_group -> real_members[j]]) + offset : (char * )(hier_group -> neighbour_sbufs[hier_group -> real_members[j]]) + offset, (do_cico) ? (char * )(vhx_module -> cico_buffer) + offset : (char * ) rbuf + offset, reduce_count, datatype);
          else // in the next levels, rnuf as eused as input sincethey contain the product of previous reductions in no CICO
            ompi_op_reduce(op, do_cico ? (char * )(vhx_module -> neighbour_cico_buffers[hier_group -> real_members[j]]) + offset : (char * )(hier_group -> neighbour_rbufs[hier_group -> real_members[j]]) + offset, (do_cico) ? (char * )(vhx_module -> cico_buffer) + offset : (char * ) rbuf + offset, reduce_count, datatype);
          reduced_chunks++;
          if (i < hier_size - 1)
            set_bytes_ready_allreduce(rank, & (vhx_module -> hier_groups[i + 1]), reduced_chunks);
        }
 
      }
	  if  ( hier_group -> real_size == 1 && i == 0){
			   if(sbuf!= MPI_IN_PLACE)
					memcpy(vhx_module->cico_buffer, sbuf, bytes_total); //TODO optimize?
			   		       set_bytes_ready_allreduce(rank, & (vhx_module -> hier_groups[i + 1]), chunks);

	  }

      if (i == hier_size - 1 && do_cico)
        vector_memcpy(rbuf, vhx_module -> cico_buffer, bytes_total,OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER); // in  CICO we need to copy the product of the final reduction to its recv buffer
      opal_atomic_wmb();
	  	//	printf("i am %d lvl %d xxx \n", rank, i);

      for (int j = 0; j < hier_group -> real_size; j++){
        hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_ack = pvt_seq;
	//	printf("i am %d lvl %d and made %d %d \n", rank, i, hier_group -> real_members[j], pvt_seq);
	  }
    } else if (hier_group -> real_members_bitmap[rank] == 1) {
		//		printf(" i am member %d lvl %d\n", rank, i);

      if (do_cico && i == 0)
        vector_memcpy(vhx_module -> cico_buffer, sbuf, bytes_total,OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
      opal_atomic_wmb();

      hier_group -> members_shared_ctrl_vars[rank].member_seq = pvt_seq;

      while (hier_group -> members_shared_ctrl_vars[rank].member_ack < pvt_seq)
		;//  printf("i am %d lvl %d pvt seq %d memeber ack %d \n", rank,  i, pvt_seq, hier_group -> members_shared_ctrl_vars[rank].member_ack);
    }

  }


}

int reduce_internal_xpmem(const void * sbuf, void * rbuf,
  int count, ompi_datatype_t * datatype, ompi_op_t * op,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module, int pvt_seq) {

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
  int hier_size = vhx_module -> hierarchy_size;
  int do_cico = (bytes_total <= OMPI_vhx_CICO_MAX);
  size_t chunks = (OMPI_vhx_CHUNK_SIZE && OMPI_vhx_CHUNK_SIZE < bytes_total) ? bytes_total / OMPI_vhx_CHUNK_SIZE : 1;
  int remainder = (OMPI_vhx_CHUNK_SIZE && OMPI_vhx_CHUNK_SIZE < bytes_total) ? bytes_total % OMPI_vhx_CHUNK_SIZE : 0;
  int has_remainder = remainder ? 1 : 0;
  init_bytes_ready_allreduce(rank, vhx_module, 0, pvt_seq);
  int complete_levels = 0;
  int first_reduction_done = 0;
  int reduced_chunks = 0;
  int sbuf_memcpy = 0;
  if (sbuf != MPI_IN_PLACE)
  ; //memcpy(rbuf, (char *)sbuf,bytes_total);
  while (reduced_chunks < chunks + has_remainder) {
    size_t offset = 0;
    int reduce_count = count;

    for (int i = 0; i < hier_size; i++) {
      vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
      hier_group -> shared_ctrl_vars[rank].rbuf_vaddr = rbuf;

      hier_group -> shared_ctrl_vars[rank].sbuf_vaddr = (void*) sbuf;
      int registration_done = 0;
      if (rank == hier_group -> leader) {
        if (i == 0)
          set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[0], chunks);

        (vhx_module -> hier_groups[i].members_shared_ctrl_vars[rank].member_seq) = pvt_seq;
   //     printf("rank %d waiting point 2 lvl %d pvt_seq is %d and member_seq is %d \n", rank, i, pvt_seq, vhx_module -> hier_groups[i].members_shared_ctrl_vars[hier_group -> real_members[hier_group -> real_size - 1]].member_seq );
        while (vhx_module -> hier_groups[i].members_shared_ctrl_vars[hier_group -> real_members[hier_group -> real_size - 1]].member_seq <pvt_seq);//{        printf("rank %d waiting point 2 lvl %d pvt_seq is %d and member_seq is %d \n", rank, i, pvt_seq, vhx_module -> hier_groups[i].members_shared_ctrl_vars[hier_group -> real_members[hier_group -> real_size - 1]].member_seq );
//}
   //     printf("rank %d waiting point 2 ended lvl %d \n", rank, i);

        int chunks_ready = 0;
        while (chunks_ready <= reduced_chunks)
          chunks_ready = __atomic_load_n( & (hier_group -> shared_ctrl_vars[hier_group -> real_members[hier_group -> real_size - 1]].reduce_ready_chunks), __ATOMIC_RELAXED);
        // 					return;
        if (i == 0 && hier_group->real_size == 1 && sbuf != MPI_IN_PLACE && sbuf_memcpy == 0){
          memcpy(rbuf, sbuf, bytes_total);
          sbuf_memcpy++;
        }
        if (i != hier_size - 1 && vhx_module -> hier_groups[i + 1].leader == rank) {
          set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[i + 1], reduced_chunks + 1);
        }
        opal_atomic_rmb();

        for (int j = 0; j < hier_group -> real_size; j++)
          hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_ack = pvt_seq;

      } else if (hier_group -> real_members_bitmap[rank] == 1) {
        int src_rank = get_left_neighbour(hier_group, rank);
      //  printf("rank %d waiting point 3 lvl %d pvt_seq is %d and member_seq is %d \n", rank, i, pvt_seq, vhx_module -> hier_groups[i].members_shared_ctrl_vars[src_rank].member_seq );

        while( (vhx_module -> hier_groups[i].members_shared_ctrl_vars[src_rank].member_seq) < pvt_seq);
		
        opal_atomic_rmb();
     //   printf("rank %d waiting point 3 ended lvl %d \n", rank, i);
   
   int chunks_ready = 0;
        do {
          chunks_ready = __atomic_load_n( & (hier_group -> shared_ctrl_vars[src_rank].reduce_ready_chunks), __ATOMIC_RELAXED);
        }
        while (chunks_ready <= reduced_chunks);

        opal_atomic_rmb();

        if (OMPI_vhx_CHUNK_SIZE) {
          offset = reduced_chunks * (OMPI_vhx_CHUNK_SIZE);

          reduce_count = (reduced_chunks == chunks && has_remainder) ? remainder : opal_min(OMPI_vhx_CHUNK_SIZE, bytes_total) / dtype_size;
        }
        //	return;
        if (i == 0 && !registration_done)
          mca_coll_vhx_get_rank_reg(src_rank, hier_group -> shared_ctrl_vars[src_rank].sbuf_vaddr,
            bytes_total, & (hier_group -> sbuf_regs[src_rank]), (mca_coll_base_module_t *)vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[src_rank]);
        else if (!registration_done)
          mca_coll_vhx_get_rank_reg(src_rank, hier_group -> shared_ctrl_vars[src_rank].rbuf_vaddr,
            bytes_total, & (hier_group -> rbuf_regs[src_rank]), (mca_coll_base_module_t *)vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[src_rank]);
        if (!registration_done)
          mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> shared_ctrl_vars[hier_group -> leader].rbuf_vaddr,
            bytes_total, & (hier_group -> rbuf_regs[hier_group -> leader]), (mca_coll_base_module_t *)vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[hier_group -> leader]);
        registration_done = 1;
        if (i == 0 && src_rank == hier_group -> leader) {


          mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> shared_ctrl_vars[hier_group -> leader].sbuf_vaddr,
            bytes_total, & (hier_group -> sbuf_regs[hier_group -> leader]), (mca_coll_base_module_t *)vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[hier_group -> leader]);
          ompi_3buff_op_reduce(op, (char * ) sbuf + offset, (char * ) hier_group -> neighbour_sbufs[hier_group -> leader] + offset, (char * ) hier_group -> neighbour_rbufs[hier_group -> leader] + offset, reduce_count, datatype);

        } else {

          ompi_op_reduce(op, (i == 0) ? (char * ) sbuf + offset : (char * ) rbuf + offset, (char * ) hier_group -> neighbour_rbufs[hier_group -> leader] + offset, reduce_count, datatype);

        }
        (vhx_module -> hier_groups[i].members_shared_ctrl_vars[rank].member_seq) = pvt_seq;

        set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[i], reduced_chunks + 1);

        //  while (hier_group -> shared_ctrl_vars[rank].coll_ack != pvt_seq);

      }
    }
    reduced_chunks++;
  }

}
int init_bytes_ready_allreduce_bcast(int my_rank, vhx_module_t * module, size_t bytes) {

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = hier_size - 1; i >= 0; i--) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
  //  	WAIT_FLAG(& hier_group -> shared_ctrl_vars[0].coll_ack, pvt_seq - 1);

    if (hier_group -> leader != my_rank)
      continue;
    volatile size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = & (hier_group -> shared_ctrl_vars[0].bytes_available);
    __atomic_store_n(tmp_bytes_available, bytes, __ATOMIC_RELAXED);

  }
  return 0;

}
int allreduce_bcast(void * buf, int count, ompi_datatype_t * datatype, int root,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module, int rank, int comm_size, int hier_size, bool do_cico) {

  vhx_module_t * vhx_module = (vhx_module_t * ) module;

//  int pvt_seq = ++(vhx_module -> pvt_coll_seq);
  size_t dtype_size, bytes_total;
  ompi_datatype_type_size(datatype, & dtype_size);
  bytes_total = count * dtype_size;

  
  if (do_cico && root == rank)
    vector_memcpy((char * )(vhx_module -> cico_buffer), buf, bytes_total,OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);

  
  if (!do_cico)
    set_vaddr(rank, module, buf);

  init_bytes_ready_allreduce_bcast(rank, vhx_module, (root == rank) ? bytes_total : 0); // the root has copied the complete buffer at the beginning of the algortihm
  int chunk_size = bytes_total;
  if (OMPI_vhx_CHUNK_SIZE)
    chunk_size = OMPI_vhx_CHUNK_SIZE;
  int src_rank;
  vhx_hier_group_t * src_hier_group = NULL;
  if (root != rank) {
    set_bcast_source(rank, (mca_coll_base_module_t * )vhx_module, & src_hier_group);
    if (!src_hier_group)
      abort();

    src_rank = src_hier_group -> leader;
  }

  
  if (root == rank) {
    

    return OMPI_SUCCESS;
  }

  size_t bytes_copied = 0;
  size_t bytes_remaining = bytes_total;
  while (bytes_remaining) {
    volatile size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = & (src_hier_group -> shared_ctrl_vars[0].bytes_available);
    size_t bytes_available = __atomic_load_n(tmp_bytes_available, __ATOMIC_RELAXED);
    opal_atomic_rmb();

    if (bytes_copied == bytes_available)
      continue;
    size_t bytes_to_be_copied = opal_min(chunk_size, bytes_remaining);

    if (do_cico)
      vector_memcpy((src_hier_group == & (vhx_module -> hier_groups[0])) ? (char * ) buf + bytes_copied : (char * )(vhx_module -> cico_buffer) + bytes_copied, (char * )(vhx_module -> neighbour_cico_buffers[src_hier_group -> leader]) + bytes_copied, bytes_to_be_copied,OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
    else {
      //	mca_coll_vhx_get_rank_reg(src_hier_group -> leader, src_hier_group -> shared_ctrl_vars[0].sbuf_vaddr,
      //	bytes_total, & (src_hier_group -> sbuf_regs[src_rank]), vhx_module, ompi_comm, & src_hier_group -> neighbour_sbufs[src_rank]);
      vector_memcpy((char * ) buf + bytes_copied, (char * )(src_hier_group -> neighbour_rbufs[src_rank]) + bytes_copied, bytes_to_be_copied,OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
    }

    opal_atomic_wmb();
    if (root != rank && (vhx_module -> hier_groups[0].leader == rank) && do_cico)
      vector_memcpy((char * ) buf + bytes_copied, (char * ) vhx_module -> cico_buffer + bytes_copied, bytes_to_be_copied,OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
    bytes_copied += bytes_to_be_copied;

    set_bytes_ready(rank, vhx_module, bytes_copied);
    opal_atomic_wmb();

    bytes_remaining -= bytes_to_be_copied;

  }


  return OMPI_SUCCESS;

}
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
  int pvt_seq = ++(vhx_module -> pvt_coll_seq);
  int hier_size = vhx_module -> hierarchy_size;
  int do_cico = (bytes_total <= OMPI_vhx_CICO_MAX);
  if (do_cico)
    reduce_internal_cico(sbuf, rbuf,
      count, datatype, op,
      ompi_comm, module, pvt_seq);
  else
    reduce_internal_xpmem(sbuf, rbuf,
     count, datatype, op,
      ompi_comm, module, pvt_seq);
  ///////BCAST/////////
  //printf("vjcskzhvks %d\n", rank);

  if (bcast) {
    allreduce_bcast(rbuf, count, datatype, 0,
      ompi_comm, module, rank, comm_size, hier_size, do_cico);

  }
  vhx_ack_wave(rank, vhx_module, pvt_seq);

  return OMPI_SUCCESS;
}

int mca_coll_vhx_allreduce(const void * sbuf, void * rbuf,
  int count, ompi_datatype_t * datatype, ompi_op_t * op,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module) {
  int rank = ompi_comm_rank(ompi_comm);

  //printf("rank %d called allreduce \n",  rank);
  return mca_coll_vhx_allreduce_internal(sbuf, rbuf,
    count, datatype, op, ompi_comm, module, true);
}
