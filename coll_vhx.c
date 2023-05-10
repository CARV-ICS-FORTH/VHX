 
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
	//	printf("ffff %p %p \n", ptr, neighbour_vaddr);
//	fflush(stdout);

	*reg = MCA_SMSC_CALL(map_peer_region, vhx_module->neighbour_endpoints[rank],
		MCA_RCACHE_FLAGS_PERSIST, neighbour_vaddr, size, &ptr);
	
	if(*reg == NULL){
		return NULL;
	}
//	printf("ff %p %d\n", ptr, *(int*)ptr);
	//fflush(stdout);
		*ptr_test= ptr;
}


int create_shmem_regions(mca_coll_base_module_t * module,
  ompi_communicator_t * comm) {

  int rank = ompi_comm_rank(comm);
  int comm_size = ompi_comm_size(comm);
  vhx_module_t * vhx_module = (vhx_module_t * ) module;

  if (rank == 0) { 
    vhx_module -> shared_ctrl_vars = vhx_shmem_create( & (vhx_module -> sync_ds), comm_size * sizeof(shared_ctrl_vars_t), "root_shared_var_ompi_vhx_component");
    for (int i = 0; i < comm_size; i++) {
      vhx_module -> shared_ctrl_vars[i].coll_seq = 0;
      vhx_module -> shared_ctrl_vars[i].coll_ack = 0;
    }
    vhx_module -> leader_cico = vhx_shmem_create( & (vhx_module -> leader_cico_ds), OMPI_vhx_CICO_MAX, "shared_cico_ompi_vhx_component");
    comm -> c_coll -> coll_barrier(comm, comm -> c_coll -> coll_barrier_module);
  } else {
    comm -> c_coll -> coll_barrier(comm, comm -> c_coll -> coll_barrier_module);
    vhx_module -> shared_ctrl_vars = vhx_shmem_attach( & (vhx_module -> sync_ds), comm_size * sizeof(shared_ctrl_vars_t), "root_shared_var_ompi_vhx_component");
    vhx_module -> leader_cico = vhx_shmem_attach( & (vhx_module -> leader_cico_ds), OMPI_vhx_CICO_MAX, "shared_cico_ompi_vhx_component");
  }

  for (int j = 0; j < vhx_module -> hierarchy_size; j++) {
    char * shmem_ctrl;

    if (rank == vhx_module -> hier_groups[j].leader) {
      opal_asprintf( & shmem_ctrl, "rank%d_level%d_shared_ctrl_var_ompi_vhx_component", rank, j);
      vhx_module -> hier_groups[j].shared_ctrl_vars = vhx_shmem_create( & (vhx_module -> hier_groups[j].sync_ds), comm_size * sizeof(shared_ctrl_vars_t), shmem_ctrl);
      free(shmem_ctrl);
      for (int i = 0; i < comm_size; i++) {
        vhx_module -> hier_groups[j].shared_ctrl_vars[i].coll_seq = 0;
        vhx_module -> hier_groups[j].shared_ctrl_vars[i].coll_ack = 0;
      }
      comm -> c_coll -> coll_barrier(comm, comm -> c_coll -> coll_barrier_module); //we block to ensure that all processes have created their shmem areas. No need to perform allgather and collect descriptors. We use same format for shmem names in all rpocs

    } else {
      comm -> c_coll -> coll_barrier(comm, comm -> c_coll -> coll_barrier_module);
      opal_asprintf( & shmem_ctrl, "rank%d_level%d_shared_ctrl_var_ompi_vhx_component", vhx_module -> hier_groups[j].leader, j);
      vhx_module -> hier_groups[j].shared_ctrl_vars = vhx_shmem_attach( & (vhx_module -> hier_groups[j].sync_ds), comm_size * sizeof(shared_ctrl_vars_t), shmem_ctrl);
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
  
  vhx_module -> initialized = true;
  return 0;
}