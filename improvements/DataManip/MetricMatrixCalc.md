# MetricMatrixCalc

## Uninitialized function pointers

### `MetricMatrixCalc.h`: calculation can call an indeterminate metric pointer

The public default constructor leaves `dp_metricFunc` uninitialized and
`calcMetricMatrix()` calls it without checking whether `setMetricFunc()` was
called. Initialize it to null and enforce the class precondition before the
loop. See
`patches/DataManip/MetricMatrixCalc/001-validate-metric-function.patch`.
