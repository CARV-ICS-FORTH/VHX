 
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

/*Currently, XPMEM support is broken since the respective code gets rewritten using differnt APIS (smsc), 
 *this legacy version is not complete and wont function properly */

void * mca_coll_vhx_get_registration(mca_rcache_base_module_t * rcache,
  void * peer_vaddr, size_t size, vhx_reg_t ** reg) {

  uint32_t flags = 0;

  uintptr_t base = (uintptr_t) peer_vaddr;

  uintptr_t bound = (uintptr_t) peer_vaddr + size;

  int ret = rcache -> rcache_register(rcache, (void * ) base, bound - base, flags,
    MCA_RCACHE_ACCESS_ANY, (mca_rcache_base_registration_t ** ) reg);

  if (ret != OPAL_SUCCESS) {
    * reg = NULL;
    return NULL;
  }

  return (void * )((uintptr_t)( * reg) -> attach_base +
    ((uintptr_t) peer_vaddr - (uintptr_t)( * reg) -> super.base));
}

static int vhx_rcache_reg_mem_cb(void * reg_data, void * base,
  size_t size, mca_rcache_base_registration_t * reg) {

  xpmem_addr_t xaddr = {
    .apid = (vhx_module_t * ) reg_data,
    .offset = (uintptr_t) base
  };
  ((vhx_reg_t * ) reg) -> attach_base = xpmem_attach(xaddr, size, NULL);

  return 0;
}

static int vhx_rcache_dereg_mem_cb(void * reg_data,
  mca_rcache_base_registration_t * reg) {

  xpmem_detach(((vhx_reg_t * ) reg) -> attach_base);
  return OPAL_SUCCESS;
}

int create_xpmem_regions(mca_coll_base_module_t * module, //Currently, XPMEME support is broken since the respective code gets rewritten using differnt APIS (smsc), this legacy version is not complete and wont dunction properly
  ompi_communicator_t * comm) {
  int rank = ompi_comm_rank(comm);
  int comm_size = ompi_comm_size(comm);
  vhx_module_t * vhx_module = (vhx_module_t * ) module;

  vhx_module -> seg_ids = malloc(sizeof(xpmem_segid_t) * comm_size);
  vhx_module -> apids = malloc(sizeof(xpmem_apid_t) * comm_size);
  xpmem_segid_t my_segid = xpmem_make(0, XPMEM_MAXADDR_SIZE,
    XPMEM_PERMIT_MODE, (void * ) 0666);

  comm -> c_coll -> coll_allgather( & (my_segid),
    sizeof(xpmem_segid_t), MPI_BYTE, & vhx_module -> seg_ids,
    sizeof(xpmem_segid_t), MPI_BYTE, comm,
    comm -> c_coll -> coll_allgather_module);
  char * rcache_name;

  opal_asprintf( & rcache_name, "vhx%d", comm_size);

  for (int i = 0; i < comm_size; i++) {
    if (i == rank) continue;

    vhx_module -> apids[i] = xpmem_get(vhx_module -> seg_ids[i],
      XPMEM_RDWR, XPMEM_PERMIT_MODE, (void * ) 0666);

    mca_rcache_base_resources_t rcache_resources = {
      .cache_name = rcache_name,
      .reg_data = & vhx_module -> apids[i],
      .sizeof_reg = sizeof(vhx_reg_t),
      .register_mem = vhx_rcache_reg_mem_cb,
      .deregister_mem = vhx_rcache_dereg_mem_cb
    };

    opal_asprintf(rcache_resources.cache_name,
      "vhx%d", i);

    vhx_module -> rcaches[i] = mca_rcache_base_module_create("grdma",
      NULL, & rcache_resources);

  }
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
int vhx_init(mca_coll_base_module_t * module,
  ompi_communicator_t * comm) {

  vhx_module_t * vhx_module = (vhx_module_t * ) module;
  vhx_data_t * data = NULL;
  int rank = ompi_comm_rank(comm);

  vhx_coll_fns_t vhx_fns = vhx_module_set_fns(comm, vhx_module -> prev_colls); //bring vanilla collective operations back for a while

  create_shmem_regions(vhx_module, comm);

  vhx_module_set_fns(comm, vhx_fns); //we replace vanilla collective operations with ours 
  //create_xpmem_regions(vhx_module, comm); // Currently unstable due to a recent code refactoring attempt TODO
  vhx_module -> initialized = true;
  return 0;
}