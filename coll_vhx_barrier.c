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
//	printf("Rank %d: level 0 coll_seq %d\n", rank, (vhx_module->hier_groups[0].shared_ctrl_vars[rank].coll_seq));
	int hier_size = vhx_module->hierarchy_size;
	//	vhx_module->hier_groups[1].shared_ctrl_vars[rank].coll_seq = pvt_seq;
	for (int i = 0; i < hier_size  ; i++){
		vhx_hier_group_t * hier_group = &(vhx_module->hier_groups[i]);
		if (i)
			hier_group->shared_ctrl_vars[rank].coll_seq = pvt_seq;
	if (rank == hier_group->leader){
		//vhx_module->shared_ctrl_vars[rank].coll_seq = ++pvt_seq;
		hier_group->shared_ctrl_vars[rank].coll_ack = pvt_seq;
		for (int j = 0; j < hier_group->real_size; j++){ //waitng for SEQ wave
		//	printf("Waiting rank %d level %d for %d to have seq %d (real size :%d\n)\n", rank, i, hier_group->real_members[j], pvt_seq, hier_group->real_size );
			while(hier_group->shared_ctrl_vars[hier_group->real_members[j]].coll_seq != pvt_seq);
			//printf("Waited rank %d level %d for %d  to have seq %d (real size %d)\n", rank, i, hier_group->real_members[j], pvt_seq, hier_group->real_size );

		}

		for (int j = 0; j <  hier_group->real_size; j++){
			hier_group->shared_ctrl_vars[hier_group->real_members[j]].coll_ack = pvt_seq;

		}

	}
	else{		

			//vhx_module->shared_ctrl_vars[rank].coll_seq = ++pvt_seq; // SEQ	
			if(hier_group->real_members_bitmap[rank]){
			//	printf("check %d pvt_seq %d\n" , rank, pvt_seq);
				while(hier_group->shared_ctrl_vars[rank].coll_ack != pvt_seq); //waiting for ACK
			}
	}
	}
	return OMPI_SUCCESS;
		}
	
