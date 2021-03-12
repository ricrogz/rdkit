/* -----------------------------------------------------------------------------
 * extend_std_map.i
 * ----------------------------------------------------------------------------- */

/**
 * SWIG 4 makes maps behave as Java HashMaps,
 * which have a put() method instead of a
 * set() one as in SWIG 3. Backport this feature.
 */
#if !defined(SWIG_VERSION) || SWIG_VERSION < 0x040000
%extend std::map {
    %rename(put) set;
}
#endif

%include "std_map.i"
