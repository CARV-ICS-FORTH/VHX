

#include "ompi_config.h"
#include "mpi.h"

#include "ompi/constants.h"
#include "ompi/datatype/ompi_datatype.h"
#include "ompi/communicator/communicator.h"
#include "opal/util/show_help.h"
#include "opal/util/minmax.h"

#include "coll_vhx.h"





int mca_coll_vhx_reduce(const void *sbuf, void *rbuf,
		int count, ompi_datatype_t *datatype, ompi_op_t *op, int root,
		ompi_communicator_t *ompi_comm, mca_coll_base_module_t *module) {
	/*reduce is treated as MPI_Allreduce without the final broadcast operation. 
	*this implies that the user *must* provide a receive buffer to the collective primitive fro ALL ranks 
	* If this requirement is not met, please use vanilla reduce instead.
	*/
	return mca_coll_vhx_allreduce_internal(sbuf, rbuf, count,
			datatype, op, ompi_comm, module, false); 
	
		}