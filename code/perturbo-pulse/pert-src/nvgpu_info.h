#include <nvml.h>

nvmlReturn_t get_gpu_uuid(unsigned int index, char *uuid, unsigned int uuid_len);

nvmlReturn_t get_gpu_meminfo(unsigned int index, unsigned long long *p_free,
    unsigned long long *p_reserved, unsigned long long *p_total,
    unsigned long long *p_used);

