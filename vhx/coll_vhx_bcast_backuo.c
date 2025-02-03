
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
int get_bcast_source(int my_rank, mca_coll_base_module_t * module){
	
	 vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i < hier_size; i++) {
	
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

  if(my_rank == hier_group->leader)
		  continue;
	  return hier_group->leader;
	  }
	  
	return -1;
	
	
}
int set_vaddr(int my_rank, mca_coll_base_module_t * module, void * sbuf) {

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i < hier_size; i++) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

    hier_group -> shared_ctrl_vars[my_rank].sbuf_vaddr = sbuf;
  }
	return 0;
}

int set_coll_seq_all_levels(int my_rank, mca_coll_base_module_t * module, int pvt_seq){ //function to change ranks's collseq in all hierarchy groups the ranks belongs to

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i < hier_size; i++) {
      vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

         (hier_group -> members_shared_ctrl_vars[my_rank].member_seq) = pvt_seq;
			if (hier_group->leader != my_rank)
				break;
  }
  return 0;

}
int init_bytes_ready(int my_rank,  vhx_module_t  * module, size_t bytes, int pvt_seq){
	
  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = hier_size -1; i >= 0; i--) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
		//	WAIT_FLAG(& hier_group -> shared_ctrl_vars[hier_group->leader].coll_ack, pvt_seq - 1);

	if (hier_group->leader != my_rank)
			continue;
    volatile  size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = &(hier_group->shared_ctrl_vars[my_rank].bytes_available);
	__atomic_store_n(tmp_bytes_available, bytes, __ATOMIC_RELAXED);

  }
  return 0;
	
}
int set_bytes_ready(int my_rank,  vhx_module_t  * module, size_t bytes){
	
  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = hier_size -1; i >= 0; i--) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
	if (hier_group->leader != my_rank)
			continue;
    volatile  size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = &(hier_group->shared_ctrl_vars[my_rank].bytes_available);
	__atomic_store_n(tmp_bytes_available, bytes, __ATOMIC_RELAXED);

  }
  return 0;
	
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

  set_leader(root, vhx_module, ompi_comm);
  if (do_cico && root == rank)
    memcpy((char * )(vhx_module -> cico_buffer), buf, bytes_total);
	
  //return OMPI_SUCCESS;
  if (!do_cico)
    set_vaddr(rank, module, buf);
						
 init_bytes_ready(rank, vhx_module, (root == rank)?bytes_total:0, pvt_seq);// the root has copied the complete buffer at the beginning of the algortihm
  for (int i = hier_size - 1; i >= 0; i--) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

    if (rank == hier_group -> leader) {
//	  opal_atomic_wmb();
     // if(OMPI_vhx_CHUNK_SIZE == 0 || hier_size == 1) // no pipelining enabled when chunk size is set to 0 by the user or when flat hierarchy is used
    //    (hier_group -> shared_ctrl_vars[rank].coll_seq) = pvt_seq;
    //  else
        set_coll_seq_all_levels(rank, vhx_module, pvt_seq );
			//printf("Level %d rank %d, \n", i, rank);

      if (i == 0 && rank != root && do_cico)  //Non root leaders need to copy data to their buf from the cico buffer in the last step
        memcpy(buf, (char * )(vhx_module -> cico_buffer), bytes_total);

    } else if (hier_group -> real_members_bitmap[rank] == 1) { //if rank belongs to hier group
		//	printf("Level %d rank %d, my leader is %d \n", i, rank, hier_group->leader);
			while(  (hier_group->members_shared_ctrl_vars[hier_group -> leader].member_seq) != pvt_seq);
		//	opal_atomic_rmb();
			if(OMPI_vhx_CHUNK_SIZE == 0 || hier_size == 1){ // no pipelining enabled when chunk size is set to 0 by the user or when flat hierarchy is used
				(hier_group -> members_shared_ctrl_vars[rank].member_seq) = pvt_seq;

				if (do_cico)
					memcpy((i == 0) ? buf : (char * )(vhx_module -> cico_buffer), (char * )(vhx_module -> neighbour_cico_buffers[hier_group -> leader]), bytes_total); // on the bottom we need to write to buf even in cico scenarios
				else {

					mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> shared_ctrl_vars[hier_group -> leader].sbuf_vaddr,
					bytes_total, & (hier_group -> sbuf_regs[hier_group -> leader]), vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[hier_group -> leader]);

					memcpy(buf, (char * )(hier_group -> neighbour_sbufs[hier_group -> leader]), bytes_total);

				}	
			//	opal_atomic_wmb();
		}
		else { //pipelining enabled

			size_t bytes_copied = 0;
			size_t bytes_remaining = bytes_total;
			while(bytes_remaining){
             volatile  size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = &(hier_group->shared_ctrl_vars[hier_group->leader].bytes_available);
				size_t bytes_available = __atomic_load_n(tmp_bytes_available, __ATOMIC_RELAXED);
                opal_atomic_rmb();

				if (bytes_copied ==  bytes_available)
					continue;
				size_t bytes_to_be_copied = opal_min(OMPI_vhx_CHUNK_SIZE, bytes_remaining);

				if (do_cico)
					memcpy((i == 0) ? (char*) buf + bytes_copied : (char * )(vhx_module -> cico_buffer) + bytes_copied, (char * )(vhx_module -> neighbour_cico_buffers[hier_group -> leader]) + bytes_copied, bytes_to_be_copied);
				else{
					mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> shared_ctrl_vars[hier_group -> leader].sbuf_vaddr,
					bytes_total, & (hier_group -> sbuf_regs[hier_group -> leader]), vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[hier_group -> leader]);
					memcpy((char*) buf + bytes_copied, (char * )(hier_group -> neighbour_sbufs[hier_group -> leader]) + bytes_copied, bytes_to_be_copied);
				}

				opal_atomic_wmb();
				bytes_copied += bytes_to_be_copied;

				set_bytes_ready(rank, vhx_module, bytes_copied);
				opal_atomic_wmb();

				set_coll_seq_all_levels(rank, vhx_module, pvt_seq );

				bytes_remaining -= bytes_to_be_copied;
							
			}
					
		}		
    }
  }

  for (int i = 0; i < hier_size; i++) {
    /*ACK Wave. Note that this ACK wave is reverse. In the barrier primitive, 
									the leader initiated the propagation of ACKs. In bcast, the group members initiate 
									while the leader waits for them*/
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
    if (rank == hier_group -> leader)
      for (int j = 0; j < hier_group -> real_size; j++) { //Waiting for the acks of the hierarchy group
        if (hier_group -> real_members[j] == rank)
          continue;
        while(  (hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_ack) != pvt_seq);
      }
    hier_group -> members_shared_ctrl_vars[rank].member_ack = pvt_seq; // we need the root's ack counter to be equal to pvt for consistency
  }

  return OMPI_SUCCESS;

}
