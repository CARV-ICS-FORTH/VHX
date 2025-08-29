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

void print_memory_decimal(const void * mem, size_t n) { //debuggubg function
	float * byte = (float * ) mem;

	for (size_t i = 0; i < n / sizeof(float); ++i) {
		printf("%f ", byte[i]);
	}

	printf("\n");
}

int vhx_get_rank_position_in_hier_group(int rank, vhx_hier_group_t * hier_group) {

	for (int i = 0; i < hier_group -> real_size; i++)
		if (hier_group -> real_members[i] == rank)
			return i;

	return -1;

}

void vhx_calculate_bytes_and_offset(int rank, vhx_hier_group_t * hier_group, int buffer_size, int data_type_size, size_t * bytes_for_rank, size_t * offset) {

	int min_bytes = 1024;
	int num_procs = hier_group -> real_size - 1;

	int eligible_procs = (buffer_size / min_bytes < num_procs) ? buffer_size / min_bytes : num_procs;

	int base_bytes = (buffer_size / eligible_procs / data_type_size) * data_type_size; //1168
	int extra_bytes = buffer_size - (base_bytes * eligible_procs); //16
	int rank_position = vhx_get_rank_position_in_hier_group(rank, hier_group) - 1;
	if( buffer_size < min_bytes){
		if(rank_position == 0){
			*bytes_for_rank = buffer_size;
			*offset = 0;
		}
		else{
			*bytes_for_rank = 0;
			*offset = 0;
		}
	return;
	}
	if (rank_position < eligible_procs) {
		* bytes_for_rank = base_bytes + (rank_position < (extra_bytes / data_type_size) ? data_type_size : 0);
	} else {
		* bytes_for_rank = 0;
	}
 
	* offset = rank_position * base_bytes + (rank_position < (extra_bytes / data_type_size) ? rank_position * data_type_size : (extra_bytes / data_type_size) * data_type_size);
	 
}
int vhx_init_bytes_ready_allreduce(int my_rank, vhx_module_t * module, size_t chunks, int pvt_seq, int comm_size) {

	vhx_module_t * vhx_module = (vhx_module_t * ) module;
	int hier_size = vhx_module -> hierarchy_size;
	for(int i = 0; i < hier_size; i++){
		vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
		hier_group->leader = hier_group->real_members[0]; //resetting default leader in case mpi_bcast has changed it
		for(int j = 0; j< comm_size; j++){
			hier_group -> neighbour_rbufs[j] = NULL;
			hier_group -> neighbour_sbufs[j] = NULL;
		}

	}
	for (int i = 0; i < hier_size; i++) {
		vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
		WAIT_FLAG( & hier_group -> shared_ctrl_vars[0].coll_ack, pvt_seq - 1);

		if (hier_group -> real_members_bitmap[my_rank] == 1) {

			volatile size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = & (hier_group -> shared_ctrl_vars[0].bytes_available);
			__atomic_store_n(tmp_bytes_available, 0, __ATOMIC_RELAXED);

		if (i == 0) {
			__atomic_store_n( & (hier_group -> shared_ctrl_vars[my_rank].reduce_ready_chunks), chunks, __ATOMIC_RELAXED);
			__atomic_store_n( & (hier_group -> members_shared_ctrl_vars[my_rank].reduce_done), 0, __ATOMIC_RELAXED);

		} else {
			__atomic_store_n( & (hier_group -> shared_ctrl_vars[my_rank].reduce_ready_chunks), 0, __ATOMIC_RELAXED);
			__atomic_store_n( & (hier_group -> members_shared_ctrl_vars[my_rank].reduce_done), 0, __ATOMIC_RELAXED);

		}
		if (hier_group -> leader != my_rank)
			break;
		}
	}
	return 0;

}

int vhx_set_bytes_ready_allreduce(int my_rank, vhx_hier_group_t * hier_group, size_t chunks) {

	__atomic_store_n( & (hier_group -> shared_ctrl_vars[my_rank].reduce_ready_chunks), chunks, __ATOMIC_RELAXED);

}

int vhx_get_left_neighbour(vhx_hier_group_t * hier_group, int rank) {
	int i = 1;
	while (1)
		if (hier_group -> real_members_bitmap[rank - i] != 1)
			i++;
		else
			break;
	return rank - i;

}
int vhx_reduce_internal_cico(const void * sbuf, void * rbuf,
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
	int first_reduction_done = -chunks -has_remainder + 1;
	vhx_init_bytes_ready_allreduce(rank, vhx_module, chunks+has_remainder, pvt_seq, comm_size);
	for (int i = 0; i < hier_size; i++) {
		vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
		hier_group -> shared_ctrl_vars[rank].rbuf_vaddr = rbuf;
		hier_group -> shared_ctrl_vars[rank].sbuf_vaddr = (char * ) sbuf;
		if (rank == hier_group -> leader) {

			hier_group -> shared_ctrl_vars[rank].coll_seq = pvt_seq;
			for (int j = 0; j < hier_group -> real_size; j++) { //waitng for SEQ wave
				if (hier_group -> real_size == 1 && i < hier_size - 1) {
					vhx_set_bytes_ready_allreduce(rank, & (vhx_module -> hier_groups[i + 1]), chunks + has_remainder);
					continue;
				}
				if (hier_group -> real_members[j] == rank)
					continue;
				WAIT_FLAG( & hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_seq, pvt_seq);
				opal_atomic_rmb();

				if (!do_cico) {
					mca_coll_vhx_get_rank_reg(hier_group -> real_members[j], hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].sbuf_vaddr,
					bytes_total, & (hier_group -> sbuf_regs[hier_group -> real_members[j]]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[hier_group -> real_members[j]]);
					mca_coll_vhx_get_rank_reg(hier_group -> real_members[j], hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].rbuf_vaddr,
					bytes_total, & (hier_group -> rbuf_regs[hier_group -> real_members[j]]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[hier_group -> real_members[j]]);
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

						reduce_count = (reduced_chunks == chunks && has_remainder) ? remainder / dtype_size : opal_min(OMPI_vhx_CHUNK_SIZE, bytes_total) / dtype_size;

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
						vhx_set_bytes_ready_allreduce(rank, & (vhx_module -> hier_groups[i + 1]), reduced_chunks);
				}

			}
			if (hier_group -> real_size == 1 && i == 0) {
        
			memcpy(vhx_module -> cico_buffer, sbuf, bytes_total); //TODO optimize?

			}

			if (i == hier_size - 1 && do_cico)
				vector_memcpy(rbuf, vhx_module -> cico_buffer, bytes_total, 1, OMPI_vhx_VECTORS_NUMBER); // in  CICO we need to copy the product of the final reduction to its recv buffer

			for (int j = 0; j < hier_group -> real_size; j++) {
				;//hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_ack = pvt_seq;
			}
		} else if (hier_group -> real_members_bitmap[rank] == 1) {

			if (do_cico && i == 0)
			vector_memcpy(vhx_module -> cico_buffer, sbuf, bytes_total, 1, OMPI_vhx_VECTORS_NUMBER);
			opal_atomic_wmb();

			hier_group -> members_shared_ctrl_vars[rank].member_seq = pvt_seq;
			}

	}

}

int vhx_check_members_for_completion(int rank, vhx_hier_group_t * hier_group, int value, int pvt_seq) {
	int group_size = hier_group -> real_size;
	int all_completed = 1;
	int is_completed = 0;
	for (int j = 0; j < group_size; j++) {
		if (rank == hier_group -> real_members[j])
			continue;
		int member_seq = hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_seq ;
		is_completed = __atomic_load_n( & (hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].reduce_done), __ATOMIC_RELAXED);

		if (is_completed < value || member_seq < pvt_seq) {
			all_completed = 0;
			break;
		}

	}
	return all_completed;
}

int vhx_reduce_internal_opt2(const void * sbuf, void * rbuf,
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
	vhx_init_bytes_ready_allreduce(rank, vhx_module, chunks+has_remainder, pvt_seq, comm_size);
	int complete_levels = 0;
	int first_reduction_done = 0;
	int reduced_chunks = 0;
	int sbuf_memcpy = 0;
  
	while (reduced_chunks < chunks + has_remainder) {
		size_t offset = 0;
		int reduce_count = count;

		for (int i = 0; i < hier_size; i++) {
			vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
			hier_group -> members_shared_ctrl_vars[rank].rbuf_vaddr = rbuf;

			hier_group -> members_shared_ctrl_vars[rank].sbuf_vaddr = (char * ) sbuf;
			int registration_done = 0;
			if (rank == hier_group -> leader) {
				if (i == 0)
					vhx_set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[0], chunks+has_remainder);

				(vhx_module -> hier_groups[i].members_shared_ctrl_vars[rank].member_seq) = pvt_seq;

				while (vhx_check_members_for_completion(rank, hier_group, reduced_chunks + 1, pvt_seq) == 0);

				if (i == 0 && hier_group -> real_size == 1  && sbuf_memcpy == 0) {
					memcpy(rbuf, sbuf, bytes_total);
					sbuf_memcpy++;
				}
				if (i != hier_size - 1) {
					vhx_set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[i + 1], reduced_chunks + 1);
				}
				opal_atomic_rmb();

				for (int j = 0; j < hier_group -> real_size; j++)
					hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_ack = pvt_seq;
			} else if (hier_group -> real_members_bitmap[rank] == 1) {
				size_t bytes_of_proc;
				size_t uniform_offset;
				if (reduced_chunks == chunks && has_remainder){
			
					vhx_calculate_bytes_and_offset(rank, hier_group, remainder, dtype_size, & bytes_of_proc, & uniform_offset);
				}
				else{
					vhx_calculate_bytes_and_offset(rank, hier_group, OMPI_vhx_CHUNK_SIZE?opal_min(bytes_total, OMPI_vhx_CHUNK_SIZE):bytes_total, dtype_size, & bytes_of_proc, & uniform_offset);
				}
				(vhx_module -> hier_groups[i].members_shared_ctrl_vars[rank].member_seq) = pvt_seq;

				opal_atomic_rmb();

				for (int j = 0; j < hier_group -> real_size; j++) {
					if (rank == hier_group -> real_members[j])
						continue;
					int src_rank = hier_group -> real_members[j];

					WAIT_FLAG( & (vhx_module -> hier_groups[i].members_shared_ctrl_vars[src_rank].member_seq), pvt_seq);
					int chunks_ready = 0;
					do {

						chunks_ready = __atomic_load_n( & (hier_group -> shared_ctrl_vars[(i==0)?hier_group -> leader:src_rank].reduce_ready_chunks), __ATOMIC_RELAXED);
					}
					while (chunks_ready <= reduced_chunks);

					opal_atomic_rmb();

					if (OMPI_vhx_CHUNK_SIZE) {
						offset = reduced_chunks * (OMPI_vhx_CHUNK_SIZE) + uniform_offset;

						reduce_count = opal_min(bytes_of_proc, bytes_total) / dtype_size;
					}

					if (i == 0 ){
						if(!hier_group -> neighbour_sbufs[src_rank])
						mca_coll_vhx_get_rank_reg(src_rank, hier_group -> members_shared_ctrl_vars[src_rank].sbuf_vaddr,
						bytes_total, & (hier_group -> sbuf_regs[src_rank]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[src_rank]);
					}
					else{
						if(!hier_group -> neighbour_rbufs[src_rank]){
							mca_coll_vhx_get_rank_reg(src_rank, hier_group -> members_shared_ctrl_vars[src_rank].rbuf_vaddr,
							bytes_total, & (hier_group -> rbuf_regs[src_rank]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[src_rank]);
				
						}
					}
				if(!hier_group -> neighbour_rbufs[hier_group -> leader])
					mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> members_shared_ctrl_vars[hier_group -> leader].rbuf_vaddr,
					bytes_total, & (hier_group -> rbuf_regs[hier_group -> leader]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[hier_group -> leader]);
				if (i == 0 && src_rank == hier_group -> leader) {
					if(!hier_group -> neighbour_sbufs[hier_group -> leader])
						mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> members_shared_ctrl_vars[hier_group -> leader].sbuf_vaddr,
						bytes_total, & (hier_group -> sbuf_regs[hier_group -> leader]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[hier_group -> leader]);
					if(reduce_count){
				
						ompi_3buff_op_reduce(op, (char * ) sbuf + offset, (char * ) hier_group -> neighbour_sbufs[hier_group -> leader] + offset, (char * ) hier_group -> neighbour_rbufs[hier_group -> leader] + offset, reduce_count, datatype);
					}
				} 
				else if (src_rank == hier_group->leader){
					if(reduce_count){
				 
					ompi_op_reduce(op, rbuf + offset, (char * ) hier_group -> neighbour_rbufs[hier_group -> leader] + offset, reduce_count, datatype);
					}
				}  
				else {

					if(reduce_count){
				 
						ompi_op_reduce(op, (i == 0) ? (char * ) hier_group -> neighbour_sbufs[src_rank] + offset : (char * ) hier_group -> neighbour_rbufs[src_rank]  + offset, (char * ) hier_group -> neighbour_rbufs[hier_group -> leader] + offset, reduce_count, datatype);
					}
				}

				}		
				opal_atomic_wmb();
				__atomic_store_n( & (hier_group -> members_shared_ctrl_vars[rank].reduce_done), reduced_chunks + 1, __ATOMIC_RELAXED);
				opal_atomic_wmb();

				vhx_set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[i], reduced_chunks + 1);

			}
		}
		reduced_chunks++;
	}
}

int vhx_reduce_internal_xpmem(const void * sbuf, void * rbuf,
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
	vhx_init_bytes_ready_allreduce(rank, vhx_module, 0, pvt_seq, comm_size);
	int complete_levels = 0;
	int first_reduction_done = 0;
	int reduced_chunks = 0;
	int sbuf_memcpy = 0;

	while (reduced_chunks < chunks + has_remainder) {
	
		size_t offset = 0;
		int reduce_count = count;

		for (int i = 0; i < hier_size; i++) {
			vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);
			hier_group -> members_shared_ctrl_vars[rank].rbuf_vaddr = rbuf;

			hier_group -> members_shared_ctrl_vars[rank].sbuf_vaddr = (char * ) sbuf;
			int registration_done = 0;
			if (rank == hier_group -> leader) {
				if (i == 0)
					vhx_set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[0], chunks+has_remainder);

				(vhx_module -> hier_groups[i].members_shared_ctrl_vars[rank].member_seq) = pvt_seq;
		
				WAIT_FLAG( & vhx_module -> hier_groups[i].members_shared_ctrl_vars[hier_group -> real_members[hier_group -> real_size - 1]].member_seq, pvt_seq);
	
				int chunks_ready = 0;
				while (chunks_ready <= reduced_chunks)
					chunks_ready = __atomic_load_n( & (hier_group -> shared_ctrl_vars[hier_group -> real_members[hier_group -> real_size - 1]].reduce_ready_chunks), __ATOMIC_RELAXED);
         
				if (i == 0 && hier_group -> real_size == 1 && sbuf_memcpy == 0) {
					memcpy(rbuf, sbuf, bytes_total);
					sbuf_memcpy++;
				}
				if (i != hier_size - 1 && vhx_module -> hier_groups[i + 1].leader == rank) {
					vhx_set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[i + 1], reduced_chunks + 1);
				}
				opal_atomic_rmb();

			for (int j = 0; j < hier_group -> real_size; j++)
				hier_group -> members_shared_ctrl_vars[hier_group -> real_members[j]].member_ack = pvt_seq;

			} else if (hier_group -> real_members_bitmap[rank] == 1) {
				int src_rank = vhx_get_left_neighbour(hier_group, rank);

				WAIT_FLAG( & (vhx_module -> hier_groups[i].members_shared_ctrl_vars[src_rank].member_seq), pvt_seq);

				opal_atomic_rmb();

				int chunks_ready = 0;
				do {
					chunks_ready = __atomic_load_n( & (hier_group -> shared_ctrl_vars[src_rank].reduce_ready_chunks), __ATOMIC_RELAXED);
				}
				while (chunks_ready <= reduced_chunks);

				opal_atomic_rmb();

				if (OMPI_vhx_CHUNK_SIZE) {
					offset = reduced_chunks * (OMPI_vhx_CHUNK_SIZE);

					reduce_count = (reduced_chunks == chunks && has_remainder) ? remainder/dtype_size : opal_min(OMPI_vhx_CHUNK_SIZE, bytes_total) / dtype_size;
				}
				if (i == 0 && !registration_done)
					mca_coll_vhx_get_rank_reg(src_rank, hier_group -> members_shared_ctrl_vars[src_rank].sbuf_vaddr,
					bytes_total, & (hier_group -> sbuf_regs[src_rank]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[src_rank]);
				else if (!registration_done)
					mca_coll_vhx_get_rank_reg(src_rank, hier_group -> members_shared_ctrl_vars[src_rank].rbuf_vaddr,
					bytes_total, & (hier_group -> rbuf_regs[src_rank]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[src_rank]);
				if (!registration_done)
					mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> members_shared_ctrl_vars[hier_group -> leader].rbuf_vaddr,
					bytes_total, & (hier_group -> rbuf_regs[hier_group -> leader]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_rbufs[hier_group -> leader]);
				registration_done = 1;
				if (i == 0 && src_rank == hier_group -> leader) {

					mca_coll_vhx_get_rank_reg(hier_group -> leader, hier_group -> members_shared_ctrl_vars[hier_group -> leader].sbuf_vaddr,
					bytes_total, & (hier_group -> sbuf_regs[hier_group -> leader]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & hier_group -> neighbour_sbufs[hier_group -> leader]);
					ompi_3buff_op_reduce(op, (char * ) sbuf + offset, (char * ) hier_group -> neighbour_sbufs[hier_group -> leader] + offset, (char * ) hier_group -> neighbour_rbufs[hier_group -> leader] + offset, reduce_count, datatype);

				} else {

					ompi_op_reduce(op, (i == 0) ? (char * ) sbuf + offset : (char * ) rbuf + offset, (char * ) hier_group -> neighbour_rbufs[hier_group -> leader] + offset, reduce_count, datatype);

				}
				(vhx_module -> hier_groups[i].members_shared_ctrl_vars[rank].member_seq) = pvt_seq;

				vhx_set_bytes_ready_allreduce(rank, & vhx_module -> hier_groups[i], reduced_chunks + 1);
		

			}
		}

		reduced_chunks++;
	}

}
int vhx_init_bytes_ready_allreduce_bcast(int my_rank, vhx_module_t * module, size_t bytes) {

	vhx_module_t * vhx_module = (vhx_module_t * ) module;
	int hier_size = vhx_module -> hierarchy_size;

	for (int i = hier_size - 1; i >= 0; i--) {
		vhx_hier_group_t * hier_group = & (vhx_module -> hier_groups[i]);

		if (hier_group -> leader != my_rank)
			continue;
		volatile size_t __attribute__((aligned(SIZEOF_SIZE_T))) * tmp_bytes_available = & (hier_group -> shared_ctrl_vars[0].bytes_available);
		__atomic_store_n(tmp_bytes_available, bytes, __ATOMIC_RELAXED);

	}
	return 0;

}
int vhx_allreduce_bcast(void * buf, int count, ompi_datatype_t * datatype, int root,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module, int rank, int comm_size, int hier_size, bool do_cico) {

	vhx_module_t * vhx_module = (vhx_module_t * ) module;

	size_t dtype_size, bytes_total;
	ompi_datatype_type_size(datatype, & dtype_size);
	bytes_total = count * dtype_size;

	if (do_cico && root == rank)
		vector_memcpy((char * )(vhx_module -> cico_buffer), buf, bytes_total, OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);

	vhx_init_bytes_ready_allreduce_bcast(rank, vhx_module, (root == rank) ? bytes_total : 0); // the root has copied the complete buffer at the beginning of the algortihm
	int chunk_size = bytes_total;
	if (OMPI_vhx_CHUNK_SIZE)
		chunk_size = OMPI_vhx_CHUNK_SIZE;
	int src_rank;
	vhx_hier_group_t * src_hier_group = NULL;
	if (root != rank) {
		vhx_set_bcast_source(rank, (mca_coll_base_module_t * ) vhx_module, & src_hier_group);
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
			vector_memcpy((src_hier_group == & (vhx_module -> hier_groups[0])) ? (char * ) buf + bytes_copied : (char * )(vhx_module -> cico_buffer) + bytes_copied, (char * )(vhx_module -> neighbour_cico_buffers[src_hier_group -> leader]) + bytes_copied, bytes_to_be_copied, OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
		else {
			mca_coll_vhx_get_rank_reg(src_hier_group -> leader, src_hier_group -> members_shared_ctrl_vars[src_hier_group -> leader].rbuf_vaddr,
			bytes_total, & (src_hier_group -> rbuf_regs[src_rank]), (mca_coll_base_module_t * ) vhx_module, ompi_comm, & src_hier_group -> neighbour_rbufs[src_rank]);

			vector_memcpy((char * ) buf + bytes_copied, (char * )(src_hier_group -> neighbour_rbufs[src_rank]) + bytes_copied, bytes_to_be_copied, OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
		}
		opal_atomic_wmb();
		if (root != rank && (vhx_module -> hier_groups[0].leader == rank) && do_cico)
			vector_memcpy((char * ) buf + bytes_copied, (char * ) vhx_module -> cico_buffer + bytes_copied, bytes_to_be_copied, OMPI_vhx_VECTOR_ELEM_SIZE, OMPI_vhx_VECTORS_NUMBER);
		bytes_copied += bytes_to_be_copied;

		vhx_set_bytes_ready(rank, vhx_module, bytes_copied);
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
	if (do_cico) {
		vhx_reduce_internal_cico(sbuf, rbuf,
		count, datatype, op,
		ompi_comm, module, pvt_seq);

	} else if (bytes_total <= OMPI_vhx_CHUNK_SIZE)
		vhx_reduce_internal_xpmem(sbuf, rbuf,
		count, datatype, op,
		ompi_comm, module, pvt_seq);
	else {
		vhx_reduce_internal_xpmem(sbuf, rbuf,
		count, datatype, op,
		ompi_comm, module, pvt_seq);

	}



	if (bcast) {
		vhx_allreduce_bcast(rbuf, count, datatype, 0,
		ompi_comm, module, rank, comm_size, hier_size, do_cico);
		vhx_ack_wave(rank, vhx_module, pvt_seq);
	} else {
		for (int i = 0; i < hier_size; i++) {

		if (vhx_module -> hier_groups[i].leader == rank)
			vhx_module -> hier_groups[i].shared_ctrl_vars -> coll_ack = pvt_seq;
		else
			WAIT_FLAG( & vhx_module -> hier_groups[i].shared_ctrl_vars -> coll_ack, pvt_seq);

		vhx_module -> hier_groups[i].members_shared_ctrl_vars[rank].member_ack = pvt_seq;
		}

	}

	return OMPI_SUCCESS;
}

int mca_coll_vhx_allreduce(const void * sbuf, void * rbuf,
  int count, ompi_datatype_t * datatype, ompi_op_t * op,
  ompi_communicator_t * ompi_comm, mca_coll_base_module_t * module) {
	int rank = ompi_comm_rank(ompi_comm);
	if(sbuf == MPI_IN_PLACE)
		sbuf = rbuf;
	return mca_coll_vhx_allreduce_internal(sbuf, rbuf,
		count, datatype, op, ompi_comm, module, true);
}
