 
 #include "coll_vhx.h"

 #define OMPI_vhx_CICO_MAX (mca_coll_vhx_component.cico_max)

 void * vhx_shmem_create(int * descriptor, size_t size,
  const char * name_chr_s) { //helper function to create shmem segmets

  char * shared_space_name = malloc(strlen(name_chr_s) + 20);
  sprintf(shared_space_name, "%s", name_chr_s);
  * descriptor = shm_open(shared_space_name, O_CREAT | O_RDWR, 0666);
  int res = ftruncate( * descriptor, size);
  free(shared_space_name);
  return mmap(0, size, PROT_WRITE | PROT_READ, MAP_SHARED, * descriptor, 0);

}

void * vhx_shmem_attach(int * descriptor, size_t size,
  const char * name_chr_s) { //helper function to attach to  shmem segmets

  char * shared_space_name = malloc(strlen(name_chr_s) + 1);
  sprintf(shared_space_name, "%s", name_chr_s);
  * descriptor = shm_open(shared_space_name, O_RDWR, 0666);
  int res = ftruncate( * descriptor, size);
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
	 abort();
	
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
      vhx_module -> neighbour_cico_buffers[i] = attach_to_cico_of_rank(i, (mca_coll_base_module_t *)vhx_module); //pre attaching to all cico bugffers of other ranks
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


// Vector memcpy function

__attribute__((target("arch=rv64gcv0p7"))) //Uncomment for vector 0.7 (or include the respective -march flags in the component's Makefile.am)
//__attribute__((target("arch=rv64gcv")))  //Uncomment for vector 1.0 (or include the respective -march flags in the component's Makefile.am)

void * vector_memcpy_e8(void * dst,
  const void * src, size_t len, unsigned long int num_vectors) {
  unsigned long int vl, iterations, remainder;
  const uint8_t * src_c = (const uint8_t * ) src;
  uint8_t * dst_c = (uint8_t * ) dst;

  if (num_vectors == 1)
    asm volatile("vsetvli %[vl], x0, e8, m1, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 2)
    asm volatile("vsetvli %[vl], x0, e8, m2, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 4)
    asm volatile("vsetvli %[vl], x0, e8, m4, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 8)
    asm volatile("vsetvli %[vl], x0, e8, m8, tu, mu": [vl]
      "=r"(vl));
  
  iterations = len / (vl);
  remainder = len % (vl);

  // Perform the vectorized copy in chunks
  for (size_t i = 0; i < iterations; ++i) {
    asm volatile(
      VLE8_V_REG("%[src]", 0) VSE8_V_REG("%[dst]", 0):
      : [src]
      "r"(src_c), [dst]
      "r"(dst_c): "memory");

    dst_c += vl;
    src_c += vl;
  }

  // Handle the remainder
  if (remainder > 0) {
    
    if (num_vectors == 1)
      asm volatile("vsetvli %[vl], %[remainder], e8, m1, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder));
    else if (num_vectors == 2)
      asm volatile("vsetvli %[vl], %[remainder], e8, m2, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder));
    else if (num_vectors == 4)
      asm volatile("vsetvli %[vl], %[remainder], e8, m4, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder));
    else if (num_vectors == 8)
      asm volatile("vsetvli %[vl], %[remainder], e8, m8, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder));

    asm volatile(
      VLE8_V_REG("%[src]", 0) VSE8_V_REG("%[dst]", 0):
      : [src]
      "r"(src_c), [dst]
      "r"(dst_c): "memory");
  }

  return dst;

}

__attribute__((target("arch=rv64gcv0p7"))) //Uncomment for vector 0.7 (or include the respective -march flags in the component's Makefile.am)
//__attribute__((target("arch=rv64gcv")))  //Uncomment for vector 1.0 (or include the respective -march flags in the component's Makefile.am)

 void * vector_memcpy_e16(void * dst,
  const void * src, size_t len, unsigned long int num_vectors) {
  unsigned long int vl, iterations, remainder;
  const uint8_t * src_c = (const uint8_t * ) src;
  uint8_t * dst_c = (uint8_t * ) dst;

  if (num_vectors == 1)
    asm volatile("vsetvli %[vl], x0, e16, m1, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 2)
    asm volatile("vsetvli %[vl], x0, e16, m2, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 4)
    asm volatile("vsetvli %[vl], x0, e16, m4, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 8)
    asm volatile("vsetvli %[vl], x0, e16, m8, tu, mu": [vl]
      "=r"(vl));
  else {
    abort();
  }
  unsigned long int vl_bytes = vl * (16 / 8);
  iterations = len / (vl_bytes);
  remainder = len % (vl_bytes);

  // Perform the vectorized copy in chunks
  for (size_t i = 0; i < iterations; ++i) {
    asm volatile(
      VLE16_V_REG("%[src]", 0) VSE16_V_REG("%[dst]", 0):
      : [src]
      "r"(src_c), [dst]
      "r"(dst_c): "memory");

    dst_c += vl_bytes;
    src_c += vl_bytes;
  }

  // Handle the remainder

  if (remainder > 0) {
    int remainder_elem = remainder / (16 / 8);
    int leftover_bytes = remainder % (16 / 8);
    if (num_vectors == 1)
      asm volatile("vsetvli %[vl], %[remainder], e16, m1, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 2)
      asm volatile("vsetvli %[vl], %[remainder], e16, m2, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 4)
      asm volatile("vsetvli %[vl], %[remainder], e16, m4, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 8)
      asm volatile("vsetvli %[vl], %[remainder], e16, m8, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    

    asm volatile(
      VLE16_V_REG("%[src]", 0) VSE16_V_REG("%[dst]", 0):
      : [src]
      "r"(src_c), [dst]
      "r"(dst_c): "memory");
    src_c += remainder_elem * 2;
    dst_c += remainder_elem * 2;
   
   if (leftover_bytes) {
      int i = leftover_bytes - 1;
      while (i >= 0) {
        dst_c[i] = src_c[i];
        i--;
      }
      return dst;

    }
  }

  return dst;

}

__attribute__((target("arch=rv64gcv0p7"))) //Uncomment for vector 0.7 (or include the respective -march flags in the component's Makefile.am)
//__attribute__((target("arch=rv64gcv")))  //Uncomment for vector 1.0 (or include the respective -march flags in the component's Makefile.am)

 void * vector_memcpy_e32(void * dst,
  const void * src, size_t len, unsigned long int num_vectors) {
  unsigned long int vl, iterations, remainder;
  const uint8_t * src_c = (const uint8_t * ) src;
  uint8_t * dst_c = (uint8_t * ) dst;

  if (num_vectors == 1)
    asm volatile("vsetvli %[vl], x0, e32, m1, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 2)
    asm volatile("vsetvli %[vl], x0, e32, m2, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 4)
    asm volatile("vsetvli %[vl], x0, e32, m4, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 8)
    asm volatile("vsetvli %[vl], x0, e32, m8, tu, mu": [vl]
      "=r"(vl));
  else {
    abort();
  }
  unsigned long int vl_bytes = vl * (32/ 8);
  iterations = len / (vl_bytes);
  remainder = len % (vl_bytes);

  // Perform the vectorized copy in chunks
  for (size_t i = 0; i < iterations; ++i) {
    asm volatile(
      VLE32_V_REG("%[src]", 0) VSE32_V_REG("%[dst]", 0):
      : [src]
      "r"(src_c), [dst]
      "r"(dst_c): "memory");

    dst_c += vl_bytes;
    src_c += vl_bytes;
  }

  // Handle the remainder

  if (remainder > 0) {
    int remainder_elem = remainder / (32 / 8);
    int leftover_bytes = remainder % (32 / 8);
    if (num_vectors == 1)
      asm volatile("vsetvli %[vl], %[remainder], e32, m1, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 2)
      asm volatile("vsetvli %[vl], %[remainder], e32, m2, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 4)
      asm volatile("vsetvli %[vl], %[remainder], e32, m4, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 8)
      asm volatile("vsetvli %[vl], %[remainder], e32, m8, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
   

    asm volatile(
      VLE32_V_REG("%[src]", 0) VSE32_V_REG("%[dst]", 0):
      : [src]
      "r"(src_c), [dst]
      "r"(dst_c): "memory");
    src_c += remainder_elem * 4;
    dst_c += remainder_elem * 4;
   
   if (leftover_bytes) {
      int i = leftover_bytes - 1;
      while (i >= 0) {
        dst_c[i] = src_c[i];
        i--;
      }
      return dst;

    }
  }

  return dst;

}

__attribute__((target("arch=rv64gcv0p7"))) //Uncomment for vector 0.7 (or include the respective -march flags in the component's Makefile.am)
//__attribute__((target("arch=rv64gcv")))  //Uncomment for vector 1.0 (or include the respective -march flags in the component's Makefile.am)

 void * vector_memcpy_e64(void * dst,
  const void * src, size_t len, unsigned long int num_vectors) {
  unsigned long int vl, iterations, remainder;
  const uint8_t * src_c = (const uint8_t * ) src;
  uint8_t * dst_c = (uint8_t * ) dst;

  if (num_vectors == 1)
    asm volatile("vsetvli %[vl], x0, e64, m1, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 2)
    asm volatile("vsetvli %[vl], x0, e64, m2, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 4)
    asm volatile("vsetvli %[vl], x0, e64, m4, tu, mu": [vl]
      "=r"(vl): "0"(vl));
  else if (num_vectors == 8)
    asm volatile("vsetvli %[vl], x0, e64, m8, tu, mu": [vl]
      "=r"(vl));
  else {
    abort();
  }
  unsigned long int vl_bytes = vl * (64/ 8);
  iterations = len / (vl_bytes);
  remainder = len % (vl_bytes);

  // Perform the vectorized copy in chunks
  for (size_t i = 0; i < iterations; ++i) {
    asm volatile(
      VLE64_V_REG("%[src]", 0) VSE64_V_REG("%[dst]", 0):
      : [src]
      "r"(src_c), [dst]
      "r"(dst_c): "memory");

    dst_c += vl_bytes;
    src_c += vl_bytes;
  }

  // Handle the remainder

  if (remainder > 0) {
    int remainder_elem = remainder / (64 / 8);
    int leftover_bytes = remainder % (64 / 8);
    if (num_vectors == 1)
      asm volatile("vsetvli %[vl], %[remainder], e64, m1, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 2)
      asm volatile("vsetvli %[vl], %[remainder], e64, m2, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 4)
      asm volatile("vsetvli %[vl], %[remainder], e64, m4, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    else if (num_vectors == 8)
      asm volatile("vsetvli %[vl], %[remainder], e64, m8, tu, mu": [vl]
        "=r"(vl): [remainder]
        "r"(remainder_elem));
    

    asm volatile(
      VLE64_V_REG("%[src]", 0) VSE64_V_REG("%[dst]", 0):
      : [src]
      "r"(src_c), [dst]
      "r"(dst_c): "memory");
    src_c += remainder_elem * 8;
    dst_c += remainder_elem * 8;
   
   if (leftover_bytes) {
      int i = leftover_bytes - 1;
      while (i >= 0) {
        dst_c[i] = src_c[i];
        i--;
      }
      return dst;

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

  create_shmem_regions((mca_coll_base_module_t *)vhx_module, comm);
  my_xpmem_init((mca_coll_base_module_t *)vhx_module, comm);
  vhx_module_set_fns(comm, vhx_fns); //we replace vanilla collective operations with ours 
  vhx_module->pvt_coll_seq = 0;
  vhx_module -> initialized = true;
  return 0;
}
