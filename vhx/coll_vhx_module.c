#include "ompi_config.h"
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "ompi/mca/coll/coll.h"
#include "ompi/mca/coll/base/base.h"
#include "opal/mca/smsc/base/base.h"

#include "opal/util/arch.h"
#include "opal/util/show_help.h"
#include "opal/util/minmax.h"

#include "coll_vhx.h"

static void mca_coll_vhx_module_construct(mca_coll_vhx_module_t *module) {

	
	module->initialized = false;
}

static void mca_coll_vhx_module_destruct(mca_coll_vhx_module_t *module) {
	if(module->initialized) {
	
	
		
		module->initialized = false;
	}
	
}

OBJ_CLASS_INSTANCE(mca_coll_vhx_module_t,
	mca_coll_base_module_t,
	mca_coll_vhx_module_construct,
	mca_coll_vhx_module_destruct);

mca_coll_base_module_t *mca_coll_vhx_module_comm_query(
		struct ompi_communicator_t *comm, int *priority) {

	if((*priority = mca_coll_vhx_component.priority) < 0)
		return NULL;
	
	if(OMPI_COMM_IS_INTER(comm) || ompi_comm_size(comm) == 1
			|| ompi_group_have_remote_peers (comm->c_local_group)) {
        
        
        return NULL;
    }
    
    int comm_size = ompi_comm_size(comm);
	for(int r = 0; r < comm_size; r++) {
		ompi_proc_t *proc = ompi_comm_peer_lookup(comm, r);
		
		if(proc->super.proc_arch != opal_local_arch) {
			
			
			return NULL;
		}
	}
	
	mca_coll_base_module_t *module =
		(mca_coll_base_module_t *) OBJ_NEW(mca_coll_vhx_module_t);
	
	if(module == NULL)
		return NULL;
	
	module->coll_module_enable = mca_coll_vhx_module_enable;
	module->coll_barrier = mca_coll_vhx_barrier;
	if(mca_smsc == NULL)
		mca_smsc_base_select();
	if(!mca_smsc_base_has_feature(MCA_SMSC_FEATURE_CAN_MAP)) {
		opal_show_help("help-coll-vhx.txt", "vhx-smsc-no-map", true);
		mca_coll_vhx_component.cico_max =  INT_MAX;
	}
	module->coll_bcast = mca_coll_vhx_bcast;
	module->coll_allreduce = mca_coll_vhx_allreduce;
	module->coll_reduce = mca_coll_vhx_reduce;
	

	
	
	
	return module;
}
#define SAVE_COLL(_dst, _f) do { \
	_dst.coll_ ## _f = comm->c_coll->coll_ ## _f; \
	_dst.coll_ ## _f ## _module = comm->c_coll->coll_ ## _f ## _module; \
	\
	if(!_dst.coll_ ## _f || !_dst.coll_ ## _f ## _module) status = OMPI_ERROR; \
	if(_dst.coll_ ## _f ## _module) OBJ_RETAIN(_dst.coll_ ## _f ## _module); \
} while(0)

#define SET_COLL(_src, _f) do { \
	comm->c_coll->coll_ ## _f = _src.coll_ ## _f; \
	comm->c_coll->coll_ ## _f ## _module = _src.coll_ ## _f ## _module; \
} while(0)
int mca_coll_vhx_module_enable(mca_coll_base_module_t *module,
		ompi_communicator_t *comm) {
	
	mca_coll_vhx_module_t *vhx_module = (mca_coll_vhx_module_t *) module;
	
	int status = OMPI_SUCCESS;
	
	SAVE_COLL(vhx_module->prev_colls, barrier);
	SAVE_COLL(vhx_module->prev_colls, bcast);
	SAVE_COLL(vhx_module->prev_colls, allreduce);
	SAVE_COLL(vhx_module->prev_colls, reduce);
	
    if(status != OMPI_SUCCESS) {
		opal_show_help("help-coll-vhx.txt", "vhx-module-save-error", true);
		return status;

	}
	if(vhx_module_prepare_hierarchy(vhx_module, comm) != 0) {
		status = OMPI_ERROR;
		return status;
		}

	if(status != OMPI_SUCCESS) {
		return status;
	}
	

	

	
	return status;
}



vhx_coll_fns_t vhx_module_set_fns(ompi_communicator_t *comm,
		vhx_coll_fns_t new) {
	
	vhx_coll_fns_t current;//////////
	int status;
	
	SAVE_COLL(current, barrier);
	SAVE_COLL(current, bcast);
	SAVE_COLL(current, allreduce);
	SAVE_COLL(current, reduce);
	
	(void) status; 
	
	SET_COLL(new, barrier);
	SET_COLL(new, bcast);
    SET_COLL(new, allreduce);

	SET_COLL(new, reduce);
	
	return current;
}


int vhx_module_prepare_hierarchy(mca_coll_vhx_module_t *module,
		ompi_communicator_t *comm) {
		int rank = ompi_comm_rank(comm);

		vhx_module_t *vhx_module = (vhx_module_t *) module;
		int comm_size = ompi_comm_size(comm);

		const char *hierarchy = mca_coll_vhx_component.hierarchy_mca;
		int hierarchy_size = 1;
		for(char * tmp =(char *)hierarchy; *tmp != '\0'; tmp++)
			if (*tmp == ',')  //counting commas in string, to determine the size of hierarchy
				hierarchy_size++;
		opal_hwloc_locality_t *  locality_levels;	
		locality_levels = malloc((hierarchy_size) * sizeof(opal_hwloc_locality_t));
		
		char * pos = strtok((char*)hierarchy, ",");
		int i=0;
		while(pos){ //we create an array of hwloc locality bit values based on the hierarchy string the user provided
			
			locality_levels[i++] = strToLoc(pos);
			pos = strtok(NULL, ",");			
		}
		if (locality_levels[i-1] != OPAL_PROC_ON_NODE){
				hierarchy_size++; //node level is not included in the user's input, however, it is implied as the top hierarchy level
				locality_levels[i] = OPAL_PROC_ON_NODE;
		}
		vhx_module->general_hierarchy = locality_levels;
		i = 0;	
		vhx_module->hier_groups = malloc( hierarchy_size * sizeof(vhx_hier_group_t));
		
		for (int j = 0; j < hierarchy_size; j++){
			int * members_bitmap = calloc(comm_size, sizeof (int));
			int * members =  malloc(comm_size *sizeof (int));
			int * real_members =  malloc(comm_size *sizeof (int));
			int * real_members_bitmap= calloc(comm_size, sizeof (int)); //in groups higher than the bottom level, only the leaders of the previous levels should be included as real members
			int k =0;
			int l = 0;
			if (j > 0){	 //in groups higher than the bottom level, only the leaders of the previous level should be included as real members
				int  member_status = ( rank == vhx_module->hier_groups[j-1].leader)?1:0;

				comm->c_coll->coll_allgather(&member_status, 1,
					MPI_INT, real_members_bitmap, 1, MPI_INT,
					comm, comm->c_coll->coll_allgather_module);
			}
			for (i =0; i<comm_size; i++){
				if (RANK_IS_LOCAL(comm, i, locality_levels[j])){
						
						members_bitmap[i] = 1;		
						members[k++] = i;

						if ( !j || (real_members_bitmap[i] == 1)){
							real_members[l++] = i;
						    real_members_bitmap[i] = 1;
						}
							
				}
			}
				vhx_module->hier_groups[j].members = members;
				vhx_module->hier_groups[j].members_bitmap = members_bitmap;
				vhx_module->hier_groups[j].real_members = real_members;
				vhx_module->hier_groups[j].real_members_bitmap = real_members_bitmap;	
				vhx_module->hier_groups[j].size = k;
				vhx_module->hier_groups[j].real_size = l;
				vhx_module->hier_groups[j].leader = real_members[0]; //member with smallest rank becomes the leader
				vhx_module->hier_groups[j].locality = locality_levels[j];


		}
		vhx_module->hierarchy_size = hierarchy_size;
		
		
		return  0;
		
		}
		
		
