#include "ompi_config.h"
#include "mpi.h"
#include "ompi/mca/mca.h"
#include "ompi/mca/coll/coll.h"
#include "ompi/constants.h"
#include "ompi/communicator/communicator.h"

#include "coll_vhx.h"


int mca_coll_vhx_barrier(ompi_communicator_t *ompi_comm,
		mca_coll_base_module_t *module){
	
	int rank = ompi_comm_rank(ompi_comm);
	int comm_size = ompi_comm_size(ompi_comm);

	vhx_module_t *vhx_module = (vhx_module_t *) module;
		if(((mca_coll_vhx_module_t *) module)->initialized == false) {
			int ret = vhx_init(module, ompi_comm);
			if(ret != 0) return OMPI_ERROR;
	}
	
	int pvt_seq = ++(vhx_module -> pvt_coll_seq);
;
	int hier_size = vhx_module->hierarchy_size;

	for (int i = 0; i < hier_size  ; i++){
		vhx_hier_group_t * hier_group = &(vhx_module->hier_groups[i]);
		
			hier_group->members_shared_ctrl_vars[rank].member_seq = pvt_seq;
	if (rank == hier_group->leader){
		hier_group->members_shared_ctrl_vars[rank].member_ack = pvt_seq;
	
	for (int j = 0; j < hier_group->real_size; j++){ //waitng for SEQ wave
		if(rank == hier_group->real_members[j] )
				continue;

			while(hier_group->members_shared_ctrl_vars[hier_group->real_members[j]].member_seq != pvt_seq);


		}
	}
	else
		break;
	}	
  for (int i = 0; i < hier_size  ; i++){
					vhx_hier_group_t * hier_group = &(vhx_module->hier_groups[i]);

	
			if((hier_group->real_members_bitmap[rank] && hier_group->leader != rank)  ){

				while(hier_group->shared_ctrl_vars->coll_ack != pvt_seq); //waiting for ACK

				
					break;
			} 
	}
	
	  for (int i = 0; i < hier_size  ; i++){
		  		vhx_hier_group_t * hier_group = &(vhx_module->hier_groups[i]);
			hier_group->members_shared_ctrl_vars[rank].member_ack = pvt_seq;
			if(hier_group->leader == rank){
				hier_group->shared_ctrl_vars->coll_ack = pvt_seq;

			}
			else
				break;
	}

	return OMPI_SUCCESS;
		}
	
