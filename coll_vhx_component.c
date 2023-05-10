/*
 * Copyright (c) 2022-2023 Computer Architecture and VLSI Systems (CARV)
 *                         Laboratory, ICS Forth. All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "ompi_config.h"
#include "mpi.h"

#include "ompi/mca/coll/coll.h"
#include "ompi/mca/coll/base/base.h"

#include "opal/mca/shmem/base/base.h"
#include "opal/util/show_help.h"

#include "coll_vhx.h"

const char *mca_coll_vhx_component_version_string =
	"Open MPI vhx collective MCA component version " OMPI_VERSION;

static int vhx_register(void);

static const char *topo_locs_str[] = {
	"node", "flat",
	"socket",
	"numa",
	"l3", "l3cache",
	"l2", "l2cache",
	"l1", "l1cache",
	"core",
	"hwthread", "thread"
};




vhx_loc_t strToLoc ( char * str){ //function to transform string to hwloc_location. Useful in the cration of hierarchy groups in vhx_module.c
	
	if (!strcasecmp(str, "node") || !strcasecmp(str, "flat"))
		return OPAL_PROC_ON_NODE;
	if (!strcasecmp(str, "socket") || !strcasecmp(str, "package"))
		return OPAL_PROC_ON_SOCKET;
	if (!strcasecmp(str, "numa"))
		return OPAL_PROC_ON_NUMA;
	if (!strcasecmp(str, "l3") || !strcasecmp(str, "l3cache"))
		return OPAL_PROC_ON_L3CACHE;
	if (!strcasecmp(str, "L2") || !strcasecmp(str, "l2cache"))
		return OPAL_PROC_ON_L2CACHE;
	if (!strcasecmp(str, "L1") || !strcasecmp(str, "l1cache"))
		return OPAL_PROC_ON_L1CACHE;
	if (!strcasecmp(str, "CORE") )
		return OPAL_PROC_ON_CORE;
	if (!strcasecmp(str, "HWTHREAD") )
		return OPAL_PROC_ON_HWTHREAD;
	 
		abort();
	
}


mca_coll_vhx_component_t mca_coll_vhx_component = {
	.super = {
		.collm_version = {
			MCA_COLL_BASE_VERSION_2_4_0,
			
			.mca_component_name = "vhx",
			MCA_BASE_MAKE_VERSION(component, OMPI_MAJOR_VERSION,
				OMPI_MINOR_VERSION, OMPI_RELEASE_VERSION),
			
			.mca_register_component_params = vhx_register,
		},  
		
		.collm_data = {
			MCA_BASE_METADATA_PARAM_CHECKPOINT
		},
		
		.collm_init_query = mca_coll_vhx_component_init_query,
		.collm_comm_query = mca_coll_vhx_module_comm_query,
	},
	
	.priority = 0,
	.print_info = false,
	.cico_max = 1024,
	.shmem_backing = NULL,

	.hierarchy_mca = "numa,socket",

	
	
	

	
};


int mca_coll_vhx_component_init_query(bool enable_progress_threads,
		bool enable_mpi_threads) {
	
	return OMPI_SUCCESS;
}



static int vhx_register(void) {
	mca_base_var_enum_t *var_enum;
	char *desc;
	int ret;
	
	/* Priority */
	
	(void) mca_base_component_var_register(&mca_coll_vhx_component.super.collm_version,
		"priority", "Priority of the vhx component",
		MCA_BASE_VAR_TYPE_INT, NULL, 0, 0, OPAL_INFO_LVL_1,
		MCA_BASE_VAR_SCOPE_READONLY, &mca_coll_vhx_component.priority);
	

	(void) mca_base_component_var_register(&mca_coll_vhx_component.super.collm_version,
		"hierarchy", desc, MCA_BASE_VAR_TYPE_STRING, NULL, 0, 0, OPAL_INFO_LVL_4,
		MCA_BASE_VAR_SCOPE_READONLY, &mca_coll_vhx_component.hierarchy_mca);
	
	free(desc); desc = NULL;
	
	
	
	mca_coll_vhx_component.shmem_backing = (access("/dev/shm", W_OK) == 0 ?
		"/dev/shm" : opal_process_info.job_session_dir);
	

	
	return OMPI_SUCCESS;
}


