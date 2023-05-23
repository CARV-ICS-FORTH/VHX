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
	
	int pvt_seq = ++(vhx_module->hier_groups[0].shared_ctrl_vars[rank].coll_seq);
	int hier_size = vhx_module->hierarchy_size;
	for (int i = 0; i < hier_size  ; i++){
		vhx_hier_group_t * hier_group = &(vhx_module->hier_groups[i]);
		if (i)
			hier_group->shared_ctrl_vars[rank].coll_seq = pvt_seq;
	if (rank == hier_group->leader){
		hier_group->shared_ctrl_vars[rank].coll_ack = pvt_seq;
	
	for (int j = 0; j < hier_group->real_size; j++){ //waitng for SEQ wave
			while(hier_group->shared_ctrl_vars[hier_group->real_members[j]].coll_seq != pvt_seq);

		}

		for (int j = 0; j <  hier_group->real_size; j++){ //The group leader starts ACK propagation
			hier_group->shared_ctrl_vars[hier_group->real_members[j]].coll_ack = pvt_seq;

		}

	}
	else{		

			if(hier_group->real_members_bitmap[rank]){
				while(hier_group->shared_ctrl_vars[rank].coll_ack != pvt_seq); //waiting for ACK
			}
	}
	}
	return OMPI_SUCCESS;
		}
	
