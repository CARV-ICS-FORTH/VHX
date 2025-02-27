
/*
 * Copyright (c) 2022-2023 Computer Architecture and VLSI Systems (CARV)
 *                         Laboratory, ICS Forth. All rights reserved.
 * $COPYRIGHT$
 *
 *
 * $HEADER$
 */

#ifndef MCA_COLL_vhx_EXPORT_H
#define MCA_COLL_vhx_EXPORT_H

#include "ompi_config.h"

#include <stdint.h>
#include <limits.h>
#include "mpi.h"
#include <string.h>
#include <fcntl.h> 
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include "ompi/mca/mca.h"
#include "ompi/mca/coll/coll.h"
#include "ompi/communicator/communicator.h"
#include "ompi/datatype/ompi_datatype.h"
#include "ompi/op/op.h"

#include "opal/mca/shmem/shmem.h"
#include "opal/mca/smsc/smsc.h"
#include "opal/mca/rcache/base/base.h"

#define RETURN_WITH_ERROR(var, err, label) do {(var) = (err); goto label;} \
	while(0)

#define OBJ_RELEASE_IF_NOT_NULL(obj) do {if((obj) != NULL) OBJ_RELEASE(obj);} while(0)

#define REALLOC(p, s, t) do {void *tmp = realloc(p, (s)*sizeof(t)); \
	if(tmp) (p) = tmp;} while(0)

#define PROC_IS_LOCAL(proc, loc) (((proc)->super.proc_flags & (loc)) == (loc))
#define RANK_IS_LOCAL(comm, rank, loc) PROC_IS_LOCAL(ompi_comm_peer_lookup(comm, rank), loc)


#define IS_SIG_ATOMIC_X_BITS(x) \
	(SIG_ATOMIC_MAX == INT ## x ## _MAX) || (SIG_ATOMIC_MAX == UINT ## x ## _MAX)

#if IS_SIG_ATOMIC_X_BITS(64)
	typedef uint64_t vhx_atomic_int;
#elif IS_SIG_ATOMIC_X_BITS(32)
	typedef uint32_t vhx_atomic_int;
#elif IS_SIG_ATOMIC_X_BITS(16)
	typedef uint16_t vhx_atomic_int;
#elif IS_SIG_ATOMIC_X_BITS(8)
	typedef uint8_t vhx_atomic_int;
#endif
// ---
typedef size_t __attribute__((aligned(SIZEOF_SIZE_T))) vhx_size_t;

// Align to cache line 
#define OMPI_vhx_CTRL_ALIGN 64

#define OMPI_vhx_OPAL_PROGRESS_CYCLE 10000

#define OMPI_vhx_CICO_MAX (mca_coll_vhx_component.cico_max)
#define OMPI_vhx_CHUNK_SIZE (mca_coll_vhx_component.chunk_size)
#define OMPI_vhx_VECTORS_NUMBER (mca_coll_vhx_component.vectors_number)
#define OMPI_vhx_VECTOR_ELEM_SIZE (mca_coll_vhx_component.vector_elem_size)


BEGIN_C_DECLS


static inline bool CHECK_FLAG(volatile int *flag,
		int thresh) {
	
	
	return  (*flag == thresh) ;
}

static inline void WAIT_FLAG(volatile int *flag,
		int thresh) {
	bool ready = false;
	
	do {
		for(int i = 0; i < 10000; i++) {
			if(CHECK_FLAG(flag, thresh)) {
				ready = true;
				break;
			}
			

		}
		
		if(!ready)
			opal_progress();
	} while(!ready);
}

// ---------------------
// ----------------------------------------

typedef opal_hwloc_locality_t vhx_loc_t;


typedef struct mca_coll_vhx_component_t mca_coll_vhx_component_t;
typedef struct mca_coll_vhx_module_t mca_coll_vhx_module_t;
typedef struct mca_coll_vhx_module_t vhx_module_t;
typedef struct xpmem_addr xpmem_addr_t;

typedef struct vhx_data_t vhx_data_t;
typedef struct vhx_comm_t vhx_comm_t;
typedef struct vhx_info vhx_info;
typedef struct vhx_comm_ctrl_t vhx_comm_ctrl_t;
typedef struct vhx_member_ctrl_t vhx_member_ctrl_t;
typedef struct shared_ctrl_vars_t shared_ctrl_vars_t;
typedef struct vhx_rank_info_t vhx_rank_info_t;
typedef struct vhx_member_info_t vhx_member_info_t;
typedef struct vhx_comm_info_t vhx_comm_info_t;

typedef void vhx_copy_data_t;
typedef int vhx_ds;
typedef struct vhx_coll_fns_t vhx_coll_fns_t;
typedef struct vhx_hier_group_t vhx_hier_group_t;
typedef void * vhx_reg_t;
OMPI_DECLSPEC extern mca_coll_vhx_component_t mca_coll_vhx_component;
OMPI_DECLSPEC OBJ_CLASS_DECLARATION(mca_coll_vhx_module_t);

// ----------------------------------------


struct vhx_coll_fns_t {
	mca_coll_base_module_allreduce_fn_t coll_allreduce;
	mca_coll_base_module_t *coll_allreduce_module;
	
	mca_coll_base_module_barrier_fn_t coll_barrier;
	mca_coll_base_module_t *coll_barrier_module;
	
	mca_coll_base_module_bcast_fn_t coll_bcast;
	mca_coll_base_module_t *coll_bcast_module;
	
	mca_coll_base_module_reduce_fn_t coll_reduce;
	mca_coll_base_module_t *coll_reduce_module;
};

struct mca_coll_vhx_component_t {
	mca_coll_base_component_t super;
	
	int priority;
	bool print_info;
	int cico_max;
	char *shmem_backing;
	char * hierarchy_mca;
	size_t chunk_size;
	int vector_copy;
	size_t vector_elem_size;
	size_t vectors_number;

};


struct mca_coll_vhx_module_t {
	mca_coll_base_module_t super;
	mca_smsc_endpoint_t ** neighbour_endpoints;
	char ** neighbour_endpoints_test;
	mca_smsc_endpoint_t * my_endpoint;
	
	bool initialized;
	vhx_ds  sync_ds;
	vhx_ds  leader_cico_ds;
	vhx_ds  cico_ds;
	vhx_ds * neighbour_cico_ds;
	int pvt_coll_seq;
	
	mca_rcache_base_module_t **rcaches;
	volatile shared_ctrl_vars_t *shared_ctrl_vars;
	vhx_hier_group_t * hier_groups;
	void *leader_cico;
	opal_hwloc_locality_t * general_hierarchy;
	int hierarchy_size;
	void * cico_buffer;
	void ** neighbour_cico_buffers;
	vhx_coll_fns_t prev_colls;
	int rank;
	
};

struct shared_ctrl_vars_t { 
	
	volatile vhx_atomic_int coll_ack __attribute__((aligned(64)));
	volatile vhx_atomic_int coll_seq __attribute__((aligned(64)));
	void* volatile sbuf_vaddr;
	void* volatile rbuf_vaddr;	
	volatile vhx_size_t bytes_available  __attribute__((aligned(64)));
	volatile vhx_size_t reduce_ready_chunks;

} __attribute__((aligned(64)));;


struct vhx_hier_group_t { 
	int leader;
	vhx_ds sync_ds; //each rank creates its own set of shared ctrl varuables. Normally only the leader of each group needs to do it but we must consider the scenario where a collective's root (chosen by the user) is not a leader
	vhx_ds members_sync_ds;
	int * members_bitmap;
	int * members;
	int * real_members;
	int * real_members_bitmap;
	int size;
	opal_hwloc_locality_t locality;
	int real_size; //equal to size if bottom level group, equal to the number of leaders in the previous level group otherwise
	int hier_level;
	vhx_ds  cico_ds;
	volatile shared_ctrl_vars_t *shared_ctrl_vars; //the set of ctrl variables for the rank
	struct vhx_member_ctrl_t *members_shared_ctrl_vars;

	void *cico_buffer;
	vhx_reg_t ** sbuf_regs ;
	void ** neighbour_sbufs;		
	vhx_reg_t **rbuf_regs;
	void ** neighbour_rbufs;

};

struct vhx_member_ctrl_t {
	volatile vhx_atomic_int member_ack; // written by member
	
	// written by member, at beginning of operation
	volatile vhx_atomic_int member_seq __attribute__((aligned(64)));
	volatile int rank;
	
	void* volatile sbuf_vaddr;
	void* volatile rbuf_vaddr;
	volatile int cico_id;
	
	// reduction progress counters, written by member
	//volatile xf_int_t reduce_ready;
	//volatile xf_int_t reduce_done;
} __attribute__((aligned(64)));

// ----------------------------------------

// coll_vhx_component.c
// --------------------

int mca_coll_vhx_component_init_query(bool enable_progress_threads,
	bool enable_mpi_threads);
vhx_loc_t strToLoc ( char * str);


// coll_vhx_module.c
// -----------------
vhx_coll_fns_t vhx_module_set_fns(ompi_communicator_t *comm,
		vhx_coll_fns_t new);
mca_coll_base_module_t *mca_coll_vhx_module_comm_query(
	ompi_communicator_t *comm, int *priority);

int mca_coll_vhx_module_enable(mca_coll_base_module_t *module,
	ompi_communicator_t *comm);


int vhx_module_prepare_hierarchy(mca_coll_vhx_module_t *module,
	ompi_communicator_t *comm);

// coll_vhx.c
// ----------

int vhx_init(mca_coll_base_module_t *module, ompi_communicator_t *comm);
void vhx_destroy_data(mca_coll_vhx_module_t *module);
void *vector_memcpy(void *dst, const void *src, size_t len, size_t element_size, unsigned long int vector_number);
void  mca_coll_vhx_get_rank_reg(int rank, void *neighbour_vaddr, size_t size, vhx_reg_t **reg, mca_coll_base_module_t * module,
  ompi_communicator_t * comm, void ** ptr_test);


// Primitives (respective file)

int mca_coll_vhx_bcast(void *buf, int count, ompi_datatype_t *datatype,
	int root, ompi_communicator_t *comm, mca_coll_base_module_t *module);

int mca_coll_vhx_barrier(ompi_communicator_t *ompi_comm,
	mca_coll_base_module_t *module);

int mca_coll_vhx_reduce(const void *sbuf, void *rbuf,
	int count, ompi_datatype_t *datatype, ompi_op_t *op, int root,
	ompi_communicator_t *comm, mca_coll_base_module_t *module);

int mca_coll_vhx_allreduce(const void *sbuf, void *rbuf,
	int count, ompi_datatype_t *datatype, ompi_op_t *op,
	ompi_communicator_t *comm, mca_coll_base_module_t *module);

int vhx_ack_wave(int rank, vhx_module_t * vhx_module, int pvt_seq);

int set_vaddr(int my_rank, mca_coll_base_module_t * module, void * sbuf);

int set_bcast_source(int my_rank, mca_coll_base_module_t * module, vhx_hier_group_t ** src_hier_group);
int set_bytes_ready(int my_rank,  vhx_module_t  * vhx_module, size_t bytes);
int mca_coll_vhx_allreduce_internal(const void * sbuf, void * rbuf,
  int count, ompi_datatype_t * datatype, ompi_op_t * op,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module, int bcast);
END_C_DECLS

#endif
