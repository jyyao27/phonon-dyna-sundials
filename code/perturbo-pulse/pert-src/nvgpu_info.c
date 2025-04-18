#include <assert.h>
#include <stdio.h>
#include <nvml.h>


/* Get the NVIDIA GPU device UUID.  See this link for more details:
 * https://docs.nvidia.com/deploy/nvml-api/group__nvmlDeviceQueries.html
 * (search for nvmlDeviceGetUUID()).
 */
nvmlReturn_t get_gpu_uuid(unsigned int index, char *uuid, unsigned int uuid_len) {
  nvmlReturn_t rc;
  nvmlDevice_t device;

  nvmlInit();

  rc = nvmlDeviceGetHandleByIndex_v2(index, &device);
  if (rc == NVML_SUCCESS)
    rc = nvmlDeviceGetUUID(device, uuid, uuid_len);

  nvmlShutdown();
  return rc;
}


/* Get the NVIDIA GPU device memory information.  See this link for more
 * details:
 * https://docs.nvidia.com/deploy/nvml-api/group__nvmlDeviceQueries.html
 * (search for nvmlDeviceGetMemoryInfo_v2()).
 */
nvmlReturn_t get_gpu_meminfo(unsigned int index, unsigned long long *p_free,
    unsigned long long *p_reserved, unsigned long long *p_total,
    unsigned long long *p_used) {
  nvmlReturn_t rc;
  nvmlDevice_t device;
  nvmlMemory_v2_t memory;

  assert(p_free != NULL);
  assert(p_reserved != NULL);
  assert(p_total != NULL);
  assert(p_used != NULL);

  nvmlInit();

  rc = nvmlDeviceGetHandleByIndex_v2(index, &device);
  if (rc != NVML_SUCCESS) {
    printf("RC1 = %d: %s\n", rc, nvmlErrorString(rc));
    goto Error;
  }

  memory.version = nvmlMemory_v2;
  rc = nvmlDeviceGetMemoryInfo_v2(device, &memory);
  if (rc != NVML_SUCCESS) {
    printf("RC2 = %d: %s\n", rc, nvmlErrorString(rc));
    goto Error;
  }

  *p_free = memory.free;
  *p_reserved = memory.reserved;
  *p_total = memory.total;
  *p_used = memory.used;

Error:

  nvmlShutdown();
  return rc;
}

