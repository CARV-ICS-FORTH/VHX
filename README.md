# VHX-OpenMPI
VHX is a developing OpenMPI component aiming to eploit XPMEM'single copy capacity, hwloc's topology awareness and Vector Instructions in order to offer optimized intra-node collective communication.


* Hierarchy

A basic feature of VHX is the construction of a n-level hierarchy based on topology information obtained through the use of hwloc and the respective API of the OPAL (Open Access Layer) library provided with OpenMPI.
 
        Hierarchy levels supported:

        - NUMA node
        - Package/Socket
        - L1/L2/L3 cache
        - Hwthread/core
        
* XPMEM/CICO support

VHX supports **zero-copy** data transfer (specifically XPMEM)
        -Implemented using OPAL's SMSC library
        -Generally used for data large message sizes (can be configured)

In addition, **copy-in-copy-out** (CICO) data transportation can be used for messages of smaller sizes or when XPMEM is not available

* Lock free synchronization

VHX does not rely on locks or atomic function to achieve synchronization between processes. Instead it utilizes a single writer-many readers scheme for the handling of control variables
  
* Data Pipeling 

VHX allows concurrent (in terms of hierarchy) transfer of data  with the use of data pipelining.

* Vector Instructions

VHX implements memcopy using Vector Instructions of the RISC-V architecture.

## Building

Before building Open MPI as described in its official documentation (<https://github.com/open-mpi/ompi/blob/master/README.md>), one must copy the component's code inside OpenMPI's ompi/mca/coll/vhx folder and execute OpenMPI's autogen.pl script (possilbe requirement of --force parameter in non developer versions of OpenMPI). Afterwards, the build process is the standard one. During OpenMPI's configuration, the user must also designate XPMEM's location with the parameter `--with-xpmem`

## Execution

When built, VHX should be a candidate betweeen other existing collective components. As a result, one must set `coll_vhx_priority` MCA param to a high value

 to get it chosen. This is achived with the following runtime paramater, which sets it to 100 (max value)
 
 `--mca coll_vhx_priority 100`
 
 Note that other candidate components must be available in the case some operations are not supported by VHX. For this reason they must be declared along with VHX using
 `--mca coll basic,libnbc,vhx`
 
 An example of an execution command follows:
   `mpirun -np 32  --mca coll basic,libnbc,vhx ----mca coll_vhx_priority 100 --mca coll_vhx_hierarchy numa --map-by core --bind-to core ./executable_name`
 
 
 
 ##  Acknowledgments

We thankfully acknowledge the support of the European Commission and the Greek
General Secretariat for Research and Innovation under the EuroHPC Programme
through the **The European-PILOT** project (GA 101034126). National contributions from the
involved state members (including the Greek General Secretariat for Research
and Innovation) match the EuroHPC funding.
