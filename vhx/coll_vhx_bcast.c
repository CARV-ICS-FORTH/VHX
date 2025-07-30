\
#include "ompi_config.h"
#include "mpi.h"
#include "ompi/constants.h"
#include "ompi/datatype/ompi_datatype.h"
#include "ompi/communicator/communicator.h"
#include "opal/util/show_help.h"
#include "opal/util/minmax.h"

#include "coll_vhx.h"

//the root becomes the leader in all hierarchy groups it participates
int vhx_set_leader(int root, mca_coll_base_module_t * module, ompi_communicator_t * ompi_comm) {

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
int vhx_set_bcast_source(int my_rank, mca_coll_base_module_t * module, vhx_hier_group_t ** src_hier_group){
	
	 vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i < hier_size; i++) {
	
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

  if(my_rank == hier_group->leader)
		  continue;
	else {
	  *src_hier_group = hier_group;
	  break;
	  }
	}
	  
	return 0;
	
	
}
int vhx_set_vaddr(int my_rank, mca_coll_base_module_t * module, void * sbuf) {

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i < hier_size; i++) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
	if(hier_group->leader == my_rank)
		hier_group -> shared_ctrl_vars[0].sbuf_vaddr = sbuf;
	else
		break;
  }
	return 0;
}

int vhx_set_coll_seq_all_levels(int my_rank, mca_coll_base_module_t * module, int pvt_seq){ //function to change ranks's collseq in all hierarchy groups the ranks belongs to

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i < hier_size; i++) {
      vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

         (hier_group -> members_shared_ctrl_vars[my_rank].member_seq) = pvt_seq;
			if (hier_group->leader != my_rank)
				break;
			else
				(hier_group -> shared_ctrl_vars[0].coll_seq) = pvt_seq;

  }
  return 0;

}
int vhx_init_bytes_ready(int my_rank,  vhx_module_t  * module, size_t bytes, int pvt_seq){
	
  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = hier_size -1; i >= 0; i--) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);


	if (hier_group->leader != my_rank)
			continue;
	 while((hier_group->shared_ctrl_vars[0].coll_ack) != pvt_seq -1);

    volatile  size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = &(hier_group->shared_ctrl_vars[0].bytes_available);
	__atomic_store_n(tmp_bytes_available, bytes, __ATOMIC_RELAXED);

  }
  return 0;
	
}
int vhx_set_bytes_ready(int my_rank,  vhx_module_t  * vhx_module, size_t bytes){
	
  int hier_size = vhx_module -> hierarchy_size;

  for (int i = 0; i <= hier_size; i++) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
	if (hier_group->leader == my_rank){
			
    volatile  size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = &(hier_group->shared_ctrl_vars[0].bytes_available);
	__atomic_store_n(tmp_bytes_available, bytes, __ATOMIC_RELAXED);
	}
	else		
		break;
  }
  return 0;
	
}

int vhx_ack_wave(int rank, vhx_module_t * vhx_module, int pvt_seq){
	
	  int hier_size = vhx_module -> hierarchy_size;
     
for(int i = 0; i < hier_size; i++) {
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
		
    hier_group -> members_shared_ctrl_vars[rank].member_ack = pvt_seq; // we need the root's ack counter to be equal to pvt for consistency

		 if (hier_group -> leader != rank)
          break;
	}
	 for (int i = 0; i < hier_size; i++) {
    /*ACK Wave. Note that this ACK wave is reverse. In the barrier primitive, 
									the leader initiated the propagation of ACKs. In bcast, the group members initiate 
									while the leader waits for them*/
    vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
    if (rank == hier_group -> leader){
      for (int j = 0; j < hier_group -> real_size; j++) { //Waiting for the acks of the hierarchy group
        if (hier_group -> real_members[j] == rank)
          continue;
	  int member_rank = hier_group -> real_members[j];
       
        while(  (hier_group -> members_shared_ctrl_vars[member_rank].member_ack) != pvt_seq);
         

      }
	      hier_group -> shared_ctrl_vars[0].coll_ack = pvt_seq; // we need the root's ack counter to be equal to pvt for consistency
	}
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

 int pvt_seq = ++(vhx_module -> pvt_coll_seq);
  size_t dtype_size, bytes_total;
  ompi_datatype_type_size(datatype, & dtype_size);
  bytes_total = count * dtype_size;

  bool do_cico = (bytes_total <= OMPI_vhx_CICO_MAX);
 
  vhx_set_leader(root, (mca_coll_base_module_t *)vhx_module, ompi_comm);
  if (do_cico && root == rank)
   vector_memcpy((char * )(vhx_module -> cico_buffer), buf, bytes_total, OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
	
  if (!do_cico)
    vhx_set_vaddr(rank, module, buf);  
						
 vhx_init_bytes_ready(rank, vhx_module, (root == rank)?bytes_total:0, pvt_seq);// the root has copied the complete buffer at the beginning of the algortihm
 int chunk_size = bytes_total;
 if(OMPI_vhx_CHUNK_SIZE)
	 chunk_size = OMPI_vhx_CHUNK_SIZE;
 int src_rank;
 vhx_hier_group_t * src_hier_group = NULL;
 if (root != rank){
	vhx_set_bcast_source(rank, (mca_coll_base_module_t *)vhx_module, &src_hier_group);
	if (!src_hier_group)
		abort();
 
  src_rank = src_hier_group->leader;
 }

			vhx_set_coll_seq_all_levels(rank, (mca_coll_base_module_t *)vhx_module, pvt_seq );
				if(root == rank){
					 vhx_ack_wave(rank, vhx_module, pvt_seq);

					 return OMPI_SUCCESS;
					 }
			WAIT_FLAG(&src_hier_group->shared_ctrl_vars[0].coll_seq, pvt_seq);
			                opal_atomic_rmb();

			size_t bytes_copied = 0;
			size_t bytes_remaining = bytes_total;
			while(bytes_remaining){
             volatile  size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = &(src_hier_group->shared_ctrl_vars[0].bytes_available);
				size_t bytes_available = __atomic_load_n(tmp_bytes_available, __ATOMIC_RELAXED);
                opal_atomic_rmb();

				if (bytes_copied ==  bytes_available)
					continue;
				size_t bytes_to_be_copied = opal_min(chunk_size, bytes_remaining);
			

				if (do_cico)
					vector_memcpy((src_hier_group == &(vhx_module->hier_groups[0])) ? (char*) buf + bytes_copied : (char * )(vhx_module -> cico_buffer) + bytes_copied, (char * )(vhx_module -> neighbour_cico_buffers[src_hier_group -> leader]) + bytes_copied, bytes_to_be_copied, OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
				else{
					mca_coll_vhx_get_rank_reg(src_hier_group -> leader, src_hier_group -> shared_ctrl_vars[0].sbuf_vaddr,
					bytes_total, & (src_hier_group -> sbuf_regs[src_rank]), (mca_coll_base_module_t *)vhx_module, ompi_comm, & src_hier_group -> neighbour_sbufs[src_rank]);
					vector_memcpy((char*) buf + bytes_copied, (char * )(src_hier_group -> neighbour_sbufs[src_rank]) + bytes_copied, bytes_to_be_copied,OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
				}

				opal_atomic_wmb();
				if(root != rank && (vhx_module->hier_groups[0].leader ==rank )&& do_cico)
					vector_memcpy((char*)buf + bytes_copied, (char*)vhx_module->cico_buffer + bytes_copied, bytes_to_be_copied,OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
					bytes_copied += bytes_to_be_copied;

				vhx_set_bytes_ready(rank, vhx_module, bytes_copied);
				opal_atomic_wmb();


				bytes_remaining -= bytes_to_be_copied;
				
			}
					
				
    
  
 vhx_ack_wave(rank, vhx_module, pvt_seq);

  return OMPI_SUCCESS;

}
