/**
This file is created by the Sandia UQTk group.
Since the inlining is turned off in dSFMT.h, the function definitions are moved from dSFMT.h to this file, leaving dSFMT.h only with declarations.
This way including dSFMT.h does not produce 'duplicate symbol' errors.
*/

#include "dSFMT.h"
#include "dsfmt_add.h"



/**
 * This function generates and returns unsigned 32-bit integer.
 * This is slower than SFMT, only for convenience usage.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  uint32_t dsfmt_genrand_uint32(dsfmt_t *dsfmt) {
    uint32_t r;
    uint64_t *psfmt64 = &dsfmt->status[0].u[0];

    if (dsfmt->idx >= DSFMT_N64) {
	dsfmt_gen_rand_all(dsfmt);
	dsfmt->idx = 0;
    }
    r = psfmt64[dsfmt->idx++] & 0xffffffffU;
    return r;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [1, 2).  This is
 * the primitive and faster than generating numbers in other ranges.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  double dsfmt_genrand_close1_open2(dsfmt_t *dsfmt) {
    double r;
    double *psfmt64 = &dsfmt->status[0].d[0];

    if (dsfmt->idx >= DSFMT_N64) {
	dsfmt_gen_rand_all(dsfmt);
	dsfmt->idx = 0;
    }
    r = psfmt64[dsfmt->idx++];
    return r;
}

/**
 * This function generates and returns unsigned 32-bit integer.
 * This is slower than SFMT, only for convenience usage.
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function.  This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  uint32_t dsfmt_gv_genrand_uint32(void) {
    return dsfmt_genrand_uint32(&dsfmt_global_data);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [1, 2).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  double dsfmt_gv_genrand_close1_open2(void) {
    return dsfmt_genrand_close1_open2(&dsfmt_global_data);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  double dsfmt_genrand_close_open(dsfmt_t *dsfmt) {
    return dsfmt_genrand_close1_open2(dsfmt) - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [0, 1).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  double dsfmt_gv_genrand_close_open(void) {
    return dsfmt_gv_genrand_close1_open2() - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1].
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function. 
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  double dsfmt_genrand_open_close(dsfmt_t *dsfmt) {
    return 2.0 - dsfmt_genrand_close1_open2(dsfmt);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1].
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  double dsfmt_gv_genrand_open_close(void) {
    return 2.0 - dsfmt_gv_genrand_close1_open2();
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  double dsfmt_genrand_open_open(dsfmt_t *dsfmt) {
    double *dsfmt64 = &dsfmt->status[0].d[0];
    union {
	double d;
	uint64_t u;
    } r;

    if (dsfmt->idx >= DSFMT_N64) {
	dsfmt_gen_rand_all(dsfmt);
	dsfmt->idx = 0;
    }
    r.d = dsfmt64[dsfmt->idx++];
    r.u |= 1;
    return r.d - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
DSFMT_PRE_INLINE  double dsfmt_gv_genrand_open_open(void) {
    return dsfmt_genrand_open_open(&dsfmt_global_data);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [1, 2) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_fill_array_close1_open2() except that this function uses
 * \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2()
 */
DSFMT_PRE_INLINE  void dsfmt_gv_fill_array_close1_open2(double array[], int size) {
    dsfmt_fill_array_close1_open2(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range (0, 1] to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() and \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
DSFMT_PRE_INLINE  void dsfmt_gv_fill_array_open_close(double array[], int size) {
    dsfmt_fill_array_open_close(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [0, 1) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
DSFMT_PRE_INLINE  void dsfmt_gv_fill_array_close_open(double array[], int size) {
    dsfmt_fill_array_close_open(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range (0, 1) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
DSFMT_PRE_INLINE  void dsfmt_gv_fill_array_open_open(double array[], int size) {
    dsfmt_fill_array_open_open(&dsfmt_global_data, array, size);
}

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 * @param dsfmt dsfmt state vector.
 * @param seed a 32-bit integer used as the seed.
 */
DSFMT_PRE_INLINE  void dsfmt_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed) {
    dsfmt_chk_init_gen_rand(dsfmt, seed, DSFMT_MEXP);
    dsfmt_reset_add();
}

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed. This function uses \b global variables.
 * @param seed a 32-bit integer used as the seed.
 * see also \sa dsfmt_init_gen_rand()
 */
DSFMT_PRE_INLINE  void dsfmt_gv_init_gen_rand(uint32_t seed) {
    dsfmt_init_gen_rand(&dsfmt_global_data, seed);
    dsfmt_reset_add();
}

/**
 * This function initializes the internal state array,
 * with an array of 32-bit integers used as the seeds.
 * @param dsfmt dsfmt state vector
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 */
DSFMT_PRE_INLINE  void dsfmt_init_by_array(dsfmt_t *dsfmt, uint32_t init_key[],
				       int key_length) {
    dsfmt_chk_init_by_array(dsfmt, init_key, key_length, DSFMT_MEXP);
}

/**
 * This function initializes the internal state array,
 * with an array of 32-bit integers used as the seeds.
 * This function uses \b global variables.
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 * see also \sa dsfmt_init_by_array()
 */
DSFMT_PRE_INLINE  void dsfmt_gv_init_by_array(uint32_t init_key[], int key_length) {
    dsfmt_init_by_array(&dsfmt_global_data, init_key, key_length);
}




#if !defined(DSFMT_DO_NOT_USE_OLD_NAMES)
/**
 * This function is just the same as dsfmt_get_idstring().
 * @return id string.
 * see also \sa dsfmt_get_idstring()
 */
DSFMT_PRE_INLINE  const char *get_idstring(void) {
    return dsfmt_get_idstring();
}

/**
 * This function is just the same as dsfmt_get_min_array_size().
 * @return minimum size of array used for fill_array functions.
 * see also \sa dsfmt_get_min_array_size()
 */
DSFMT_PRE_INLINE  int get_min_array_size(void) {
    return dsfmt_get_min_array_size();
}

/**
 * This function is just the same as dsfmt_gv_init_gen_rand().
 * @param seed a 32-bit integer used as the seed.
 * see also \sa dsfmt_gv_init_gen_rand(), \sa dsfmt_init_gen_rand().
 */
DSFMT_PRE_INLINE  void init_gen_rand(uint32_t seed) {
    dsfmt_gv_init_gen_rand(seed);
    dsfmt_reset_add();
}

/**
 * This function is just the same as dsfmt_gv_init_by_array().
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 * see also \sa dsfmt_gv_init_by_array(), \sa dsfmt_init_by_array().
 */
DSFMT_PRE_INLINE  void init_by_array(uint32_t init_key[], int key_length) {
    dsfmt_gv_init_by_array(init_key, key_length);
}

/**
 * This function is just the same as dsfmt_gv_genrand_close1_open2().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_close1_open2() \sa
 * dsfmt_gv_genrand_close1_open2()
 */
DSFMT_PRE_INLINE  double genrand_close1_open2(void) {
    return dsfmt_gv_genrand_close1_open2();
}

/**
 * This function is just the same as dsfmt_gv_genrand_close_open().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_close_open() \sa
 * dsfmt_gv_genrand_close_open()
 */
DSFMT_PRE_INLINE  double genrand_close_open(void) {
    return dsfmt_gv_genrand_close_open();
}

/**
 * This function is just the same as dsfmt_gv_genrand_open_close().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_open_close() \sa
 * dsfmt_gv_genrand_open_close()
 */
DSFMT_PRE_INLINE  double genrand_open_close(void) {
    return dsfmt_gv_genrand_open_close();
}

/**
 * This function is just the same as dsfmt_gv_genrand_open_open().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_open_open() \sa
 * dsfmt_gv_genrand_open_open()
 */
DSFMT_PRE_INLINE  double genrand_open_open(void) {
    return dsfmt_gv_genrand_open_open();
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_open_close().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_open_close(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
DSFMT_PRE_INLINE  void fill_array_open_close(double array[], int size) {
    dsfmt_gv_fill_array_open_close(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_close_open().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_close_open(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
DSFMT_PRE_INLINE void fill_array_close_open(double array[], int size) {
    dsfmt_gv_fill_array_close_open(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_open_open().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_open_open(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
DSFMT_PRE_INLINE void fill_array_open_open(double array[], int size) {
    dsfmt_gv_fill_array_open_open(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_close1_open2().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
DSFMT_PRE_INLINE  void fill_array_close1_open2(double array[], int size) {
    dsfmt_gv_fill_array_close1_open2(array, size);
}
#endif /* DSFMT_DO_NOT_USE_OLD_NAMES */
