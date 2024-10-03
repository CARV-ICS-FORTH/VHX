 
 #include "coll_vhx.h"

 #define OMPI_vhx_CICO_MAX (mca_coll_vhx_component.cico_max)

 void * vhx_shmem_create(int * descriptor, size_t size,
  const char * name_chr_s) { //helper function to create shmem segmets

  char * shared_space_name = malloc(strlen(name_chr_s) + 20);
  sprintf(shared_space_name, "%s", name_chr_s);
  * descriptor = shm_open(shared_space_name, O_CREAT | O_RDWR, 0666);
  ftruncate( * descriptor, size);
  free(shared_space_name);
  return mmap(0, size, PROT_WRITE | PROT_READ, MAP_SHARED, * descriptor, 0);

}

void * vhx_shmem_attach(int * descriptor, size_t size,
  const char * name_chr_s) { //helper function to attach to  shmem segmets

  char * shared_space_name = malloc(strlen(name_chr_s) + 1);
  sprintf(shared_space_name, "%s", name_chr_s);
  * descriptor = shm_open(shared_space_name, O_RDWR, 0666);
  ftruncate( * descriptor, size);
  free(shared_space_name);
  return mmap(0, size, PROT_WRITE | PROT_READ, MAP_SHARED, * descriptor, 0);

}

void * attach_to_cico_of_rank(int rank, mca_coll_base_module_t * module) { ///helper function to attach to cico uffer of another rank

  char * shmem_cico;
  vhx_module_t * vhx_module = (vhx_module_t * ) module;

  opal_asprintf( & shmem_cico, "rank%d_shared_cico_ompi_vhx_component", rank);
  return vhx_shmem_attach( & (vhx_module -> neighbour_cico_ds[rank]), OMPI_vhx_CICO_MAX, shmem_cico);

} 


void  mca_coll_vhx_get_rank_reg(int rank, void *neighbour_vaddr, size_t size, vhx_reg_t **reg, mca_coll_base_module_t * module, 
  ompi_communicator_t * comm, void ** ptr_test) {

	vhx_module_t * vhx_module = (vhx_module_t * ) module;
		
	ompi_proc_t *proc =  ompi_comm_peer_lookup(comm, rank);
	 if(vhx_module->neighbour_endpoints[rank] == NULL) {
		
		vhx_module->neighbour_endpoints[rank] = MCA_SMSC_CALL(get_endpoint, &(proc->super));
	 }
	

	if(vhx_module->neighbour_endpoints[rank] == NULL)
	 return NULL;
	
	void *ptr;


	*reg = MCA_SMSC_CALL(map_peer_region, vhx_module->neighbour_endpoints[rank],
		MCA_RCACHE_FLAGS_PERSIST, neighbour_vaddr, size, &ptr);
	if(*reg == NULL){
		abort();;
	}

		*ptr_test= ptr;
}


int create_shmem_regions(mca_coll_base_module_t * module,
  ompi_communicator_t * comm) {

  int rank = ompi_comm_rank(comm);
  int comm_size = ompi_comm_size(comm);
  vhx_module_t * vhx_module = (vhx_module_t * ) module;


  for (int j = 0; j < vhx_module -> hierarchy_size; j++) {
    char * shmem_ctrl;

    if (rank == vhx_module -> hier_groups[j].leader) {
      opal_asprintf( & shmem_ctrl, "rank%d_level%d_shared_ctrl_var_ompi_vhx_component", rank, j);
      vhx_module -> hier_groups[j].shared_ctrl_vars = vhx_shmem_create( & (vhx_module -> hier_groups[j].sync_ds), comm_size * sizeof(shared_ctrl_vars_t), shmem_ctrl);
	  printf("shared_ctrl vars initialiazed size :%d \n", comm_size );
	  printf("rank %d shared ctrl var 0: %p 1: %p 2: %p 3: %p \n", rank, &(vhx_module -> hier_groups[j].shared_ctrl_vars[0]), &(vhx_module -> hier_groups[j].shared_ctrl_vars[1]), &(vhx_module -> hier_groups[j].shared_ctrl_vars[2]),&( vhx_module -> hier_groups[j].shared_ctrl_vars[3]));
      free(shmem_ctrl);
	  opal_asprintf( & shmem_ctrl, "rank%d_level%d_members_shared_ctrl_var_ompi_vhx_component", rank, j);
      vhx_module -> hier_groups[j].members_shared_ctrl_vars = vhx_shmem_create( & (vhx_module -> hier_groups[j].members_sync_ds), comm_size * sizeof(vhx_member_ctrl_t), shmem_ctrl);
      free(shmem_ctrl);
      for (int i = 0; i < comm_size; i++) {
        vhx_module -> hier_groups[j].shared_ctrl_vars[i].coll_seq = 0;
        vhx_module -> hier_groups[j].shared_ctrl_vars[i].coll_ack = 0;
		vhx_module -> hier_groups[j].members_shared_ctrl_vars[i].member_seq = 0; 
        vhx_module -> hier_groups[j].members_shared_ctrl_vars[i].member_ack = 0;
      }
      comm -> c_coll -> coll_barrier(comm, comm -> c_coll -> coll_barrier_module); //we block to ensure that all processes have created their shmem areas. No need to perform allgather and collect descriptors. We use same format for shmem names in all rpocs

    } else {
      comm -> c_coll -> coll_barrier(comm, comm -> c_coll -> coll_barrier_module);
      opal_asprintf( & shmem_ctrl, "rank%d_level%d_shared_ctrl_var_ompi_vhx_component", vhx_module -> hier_groups[j].leader, j);
      vhx_module -> hier_groups[j].shared_ctrl_vars = vhx_shmem_attach( & (vhx_module -> hier_groups[j].sync_ds), comm_size * sizeof(shared_ctrl_vars_t), shmem_ctrl);
      free(shmem_ctrl);
	  opal_asprintf( & shmem_ctrl, "rank%d_level%d_members_shared_ctrl_var_ompi_vhx_component", vhx_module -> hier_groups[j].leader, j);
      vhx_module -> hier_groups[j].members_shared_ctrl_vars = vhx_shmem_attach( & (vhx_module -> hier_groups[j].members_sync_ds), comm_size * sizeof(vhx_member_ctrl_t), shmem_ctrl);
      free(shmem_ctrl);
    }
    char * shmem_cico;

    opal_asprintf( & shmem_cico, "rank%d_shared_cico_ompi_vhx_component", rank);
    vhx_module -> cico_buffer = vhx_shmem_create( & (vhx_module -> cico_ds), OMPI_vhx_CICO_MAX, shmem_cico);
    vhx_module -> neighbour_cico_buffers = malloc(sizeof(char * ) * comm_size);
    vhx_module -> neighbour_cico_ds = malloc(sizeof(vhx_ds) * comm_size);

    free(shmem_cico);
    comm -> c_coll -> coll_barrier(comm, comm -> c_coll -> coll_barrier_module);

    for (int i = 0; i < comm_size; i++) {
      if (i == rank)
        continue;
      vhx_module -> neighbour_cico_buffers[i] = attach_to_cico_of_rank(i, vhx_module); //pre attaching to all cico bugffers of other ranks
    }

  }
  return 0;
}

int my_xpmem_init(mca_coll_base_module_t * module,
  ompi_communicator_t * comm){
	vhx_module_t * vhx_module = (vhx_module_t * ) module;

  int comm_size = ompi_comm_size(comm);
  vhx_module->neighbour_endpoints = malloc(sizeof(char * ) * comm_size);
  if (vhx_module->neighbour_endpoints == NULL)
	 abort();
  for (int i = 0; i < comm_size; i++){
	vhx_module->neighbour_endpoints[i] = NULL;

  }
	 for (int j = 0; j < vhx_module -> hierarchy_size; j++) {
		vhx_module->hier_groups[j].neighbour_rbufs = malloc (sizeof(void *) * comm_size);
		vhx_module->hier_groups[j].neighbour_sbufs = malloc (sizeof(void *) * comm_size);
		vhx_module->hier_groups[j].sbuf_regs =  malloc (sizeof(vhx_reg_t *) * comm_size);
		vhx_module->hier_groups[j].rbuf_regs =  malloc (sizeof(vhx_reg_t *) * comm_size);

	 }
}


#if defined(__riscv) 

#define xstr(s) str(s)
#define str(s) #s


#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// Macro to generate vector load and store instructions
#define VLE8_V_REG(src, reg) "vle8.v v" STR(reg) ", (" src ")\n\t"
#define VSE8_V_REG(dst, reg) "vse8.v v" STR(reg) ", (" dst ")\n\t"

#define VLE16_V_REG(src, reg) "vle16.v v" STR(reg) ", (" src ")\n\t"
#define VSE16_V_REG(dst, reg) "vse16.v v" STR(reg) ", (" dst ")\n\t"

#define VLE32_V_REG(src, reg) "vle32.v v" STR(reg) ", (" src ")\n\t"
#define VSE32_V_REG(dst, reg) "vse32.v v" STR(reg) ", (" dst ")\n\t"

#define VLE64_V_REG(src, reg) "vle64.v v" STR(reg) ", (" src ")\n\t"
#define VSE64_V_REG(dst, reg) "vse64.v v" STR(reg) ", (" dst ")\n\t"

#define VLE128_V_REG(src, reg) "vle128.v v" STR(reg) ", (" src ")\n\t"
#define VSE128_V_REG(dst, reg) "vse128.v v" STR(reg) ", (" dst ")\n\t"


// Vector memcpy function


void *vector_memcpy_e8(void *dst, const void *src, size_t len, unsigned long int num_vectors){
    unsigned long int vl, iterations, remainder;
    const uint8_t *src_c = (const uint8_t *)src;
    uint8_t *dst_c = (uint8_t *)dst;
  

   
    asm volatile ("vsetvli %[vl], zero, e8, m1, ta, ma"
                  : [vl] "=r" (vl)
                  : "0" (vl));
    
    
    iterations = len / (vl * num_vectors);
    remainder = len % (vl * num_vectors);

    // Perform the vectorized copy in chunks
    for (size_t i = 0; i < iterations; ++i) {
        for (int j = 0; j < num_vectors; ++j) {
            switch (j) {
                case 0:
                    asm volatile (
                        VLE8_V_REG("%[src]", 2)
                        VSE8_V_REG("%[dst]", 2)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 1:
                    asm volatile (
                        VLE8_V_REG("%[src]", 3)
                        VSE8_V_REG("%[dst]", 3)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 2:
                    asm volatile (
                        VLE8_V_REG("%[src]", 4)
                        VSE8_V_REG("%[dst]", 4)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 3:
                    asm volatile (
                        VLE8_V_REG("%[src]", 5)
                        VSE8_V_REG("%[dst]", 5)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
            }
        }
        dst_c += vl * num_vectors;
        src_c += vl * num_vectors;
    }

    // Handle the remainder
    if (remainder > 0) {
        asm volatile ("vsetvli %[vl], %[remainder], e8, m1, ta, ma"
                      : [vl] "=r" (vl)
                      : [remainder] "r" (remainder));
        for (int j = 0; j < num_vectors && j * vl < remainder; ++j) {
            switch (j) {
                case 0:
                    asm volatile (
                        VLE8_V_REG("%[src]", 2)
                        VSE8_V_REG("%[dst]", 2)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 1:
                    asm volatile (
                        VLE8_V_REG("%[src]", 3)
                        VSE8_V_REG("%[dst]", 3)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 2:
                    asm volatile (
                        VLE8_V_REG("%[src]", 4)
                        VSE8_V_REG("%[dst]", 4)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 3:
                    asm volatile (
                        VLE8_V_REG("%[src]", 5)
                        VSE8_V_REG("%[dst]", 5)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
            }
        }
    }

    return dst;
}


void *vector_memcpy_e16(void *dst, const void *src, size_t len, unsigned long int num_vectors){
  unsigned long int vl, iterations, remainder;
    const uint8_t *src_c = (const uint8_t *)src;
    uint8_t *dst_c = (uint8_t *)dst;
  
    asm volatile ("vsetvli %[vl], zero, e16, m1, ta, ma"
                  : [vl] "=r" (vl)
                  : "0" (vl));

    
    iterations = len / (vl * num_vectors);
    remainder = len % (vl * num_vectors);

    // Perform the vectorized copy in chunks
    for (size_t i = 0; i < iterations; ++i) {
        for (int j = 0; j < num_vectors; ++j) {
            switch (j) {
                case 0:
                    asm volatile (
                        VLE16_V_REG("%[src]", 2)
                        VSE16_V_REG("%[dst]", 2)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 1:
                    asm volatile (
                        VLE16_V_REG("%[src]", 3)
                        VSE16_V_REG("%[dst]", 3)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 2:
                    asm volatile (
                        VLE16_V_REG("%[src]", 4)
                        VSE16_V_REG("%[dst]", 4)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 3:
                    asm volatile (
                        VLE16_V_REG("%[src]", 5)
                        VSE16_V_REG("%[dst]", 5)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
            }
        }
        dst_c += vl * num_vectors;
        src_c += vl * num_vectors;
    }

    // Handle the remainder
    if (remainder > 0) {
        asm volatile ("vsetvli %[vl], %[remainder], e8, m1, ta, ma"
                      : [vl] "=r" (vl)
                      : [remainder] "r" (remainder));
        for (int j = 0; j < num_vectors && j * vl < remainder; ++j) {
            switch (j) {
                case 0:
                    asm volatile (
                        VLE16_V_REG("%[src]", 2)
                        VSE16_V_REG("%[dst]", 2)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 1:
                    asm volatile (
                        VLE16_V_REG("%[src]", 3)
                        VSE16_V_REG("%[dst]", 3)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 2:
                    asm volatile (
                        VLE16_V_REG("%[src]", 4)
                        VSE16_V_REG("%[dst]", 4)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 3:
                    asm volatile (
                        VLE16_V_REG("%[src]", 5)
                        VSE16_V_REG("%[dst]", 5)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
            }
        }
    }

    return dst;
}

void *vector_memcpy_e32(void *dst, const void *src, size_t len, unsigned long int num_vectors){

	  unsigned long int vl, iterations, remainder;
    const uint8_t *src_c = (const uint8_t *)src;
    uint8_t *dst_c = (uint8_t *)dst;
  
    asm volatile ("vsetvli %[vl], zero, e32, m1, ta, ma"
                  : [vl] "=r" (vl)
                  : "0" (vl));

    // Calculate iterations and remainder
    iterations = len / (vl * num_vectors);
    remainder = len % (vl * num_vectors);

    // Perform the vectorized copy in chunks
    for (size_t i = 0; i < iterations; ++i) {
        for (int j = 0; j < num_vectors; ++j) {
            switch (j) {
                case 0:
                    asm volatile (
                        VLE32_V_REG("%[src]", 2)
                        VSE32_V_REG("%[dst]", 2)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 1:
                    asm volatile (
                        VLE32_V_REG("%[src]", 3)
                        VSE32_V_REG("%[dst]", 3)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 2:
                    asm volatile (
                        VLE32_V_REG("%[src]", 4)
                        VSE32_V_REG("%[dst]", 4)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 3:
                    asm volatile (
                        VLE32_V_REG("%[src]", 5)
                        VSE32_V_REG("%[dst]", 5)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
            }
        }
        dst_c += vl * num_vectors;
        src_c += vl * num_vectors;
    }

    // Handle the remainder
    if (remainder > 0) {
        asm volatile ("vsetvli %[vl], %[remainder], e8, m1, ta, ma"
                      : [vl] "=r" (vl)
                      : [remainder] "r" (remainder));
        for (int j = 0; j < num_vectors && j * vl < remainder; ++j) {
            switch (j) {
                case 0:
                    asm volatile (
                        VLE32_V_REG("%[src]", 2)
                        VSE32_V_REG("%[dst]", 2)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 1:
                    asm volatile (
                        VLE32_V_REG("%[src]", 3)
                        VSE32_V_REG("%[dst]", 3)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 2:
                    asm volatile (
                        VLE32_V_REG("%[src]", 4)
                        VSE32_V_REG("%[dst]", 4)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 3:
                    asm volatile (
                        VLE32_V_REG("%[src]", 5)
                        VSE32_V_REG("%[dst]", 5)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
            }
        }
    }

    return dst;

}
void *vector_memcpy_e64(void *dst, const void *src, size_t len, unsigned long int num_vectors){

	  unsigned long int vl, iterations, remainder;
    const uint8_t *src_c = (const uint8_t *)src;
    uint8_t *dst_c = (uint8_t *)dst;
  
    asm volatile ("vsetvli %[vl], zero, e64, m1, ta, ma"
                  : [vl] "=r" (vl)
                  : "0" (vl));

    // Calculate iterations and remainder
    iterations = len / (vl * num_vectors);
    remainder = len % (vl * num_vectors);

    // Perform the vectorized copy in chunks
    for (size_t i = 0; i < iterations; ++i) {
        for (int j = 0; j < num_vectors; ++j) {
            switch (j) {
                case 0:
                    asm volatile (
                        VLE64_V_REG("%[src]", 2)
                        VSE64_V_REG("%[dst]", 2)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 1:
                    asm volatile (
                        VLE64_V_REG("%[src]", 3)
                        VSE64_V_REG("%[dst]", 3)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 2:
                    asm volatile (
                        VLE64_V_REG("%[src]", 4)
                        VSE64_V_REG("%[dst]", 4)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 3:
                    asm volatile (
                        VLE64_V_REG("%[src]", 5)
                        VSE64_V_REG("%[dst]", 5)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
            }
        }
        dst_c += vl * num_vectors;
        src_c += vl * num_vectors;
    }

    // Handle the remainder
    if (remainder > 0) {
        asm volatile ("vsetvli %[vl], %[remainder], e8, m1, ta, ma"
                      : [vl] "=r" (vl)
                      : [remainder] "r" (remainder));
        for (int j = 0; j < num_vectors && j * vl < remainder; ++j) {
            switch (j) {
                case 0:
                    asm volatile (
                        VLE64_V_REG("%[src]", 2)
                        VSE64_V_REG("%[dst]", 2)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 1:
                    asm volatile (
                        VLE64_V_REG("%[src]", 3)
                        VSE64_V_REG("%[dst]", 3)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 2:
                    asm volatile (
                        VLE64_V_REG("%[src]", 4)
                        VSE64_V_REG("%[dst]", 4)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
                case 3:
                    asm volatile (
                        VLE64_V_REG("%[src]", 5)
                        VSE64_V_REG("%[dst]", 5)
                        :
                        : [src] "r" (src_c + j * vl), [dst] "r" (dst_c + j * vl)
                        : "memory");
                    break;
            }
        }
    }

    return dst;

}

 
void *vector_memcpy(void *dst, const void *src, size_t len, size_t element_size, unsigned long int vector_number) {

        if(element_size == 1)
                return vector_memcpy_e8(dst, src,  len,   vector_number);
        else if (element_size == 2)
                return vector_memcpy_e16(dst, src,  len, vector_number);
        else if (element_size == 4)
                return vector_memcpy_e32(dst, src,  len, vector_number);
        else if (element_size == 8)
                return vector_memcpy_e64(dst, src,  len, vector_number);
      
	else if (element_size == 0 || vector_number == 0)
		return memcpy(dst, src, len);


    //printf("dddddd\n");
    //fflush(stdout);

}
#else
// Fallback to standard memcpy
#include <string.h>
void *vector_memcpy(void *dst, const void *src, size_t len, size_t element_size, unsigned long int vector_number) {
    return memcpy(dst, src, len);
}

#endif


int vhx_init(mca_coll_base_module_t * module,
  ompi_communicator_t * comm) {

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  vhx_data_t * data = NULL;
  int rank = ompi_comm_rank(comm);
  int comm_size = ompi_comm_size(comm);

  vhx_coll_fns_t vhx_fns = vhx_module_set_fns(comm, vhx_module -> prev_colls); //bring vanilla collective operations back for a while

  create_shmem_regions(vhx_module, comm);
  my_xpmem_init(vhx_module, comm);
  vhx_module_set_fns(comm, vhx_fns); //we replace vanilla collective operations with ours 
  vhx_module->pvt_coll_seq = 0;
  vhx_module -> initialized = true;
  return 0;
}
