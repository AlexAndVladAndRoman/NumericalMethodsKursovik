import java.util.Random

/**
 * System is represented by an array of functions that returns 0 on the answer.
 */


object RandomHolder {
    val random = Random(System.currentTimeMillis())
}

class EquationSystem(private val functions: Array<Function>) {
    private val n: Int = functions.size

    /**
     * @param x argument
     * @return discrepancy = sum fi(x)*fi(x), i = 0..n-1
     */
    fun discrepancy(x: DoubleArray): Double {
        var discrepancy = 0.0
        var q: Double
        for (i in 0 until n) {
            q = functions[i].calculate(x)
            discrepancy += q * q
        }
        return discrepancy
    }

    /**
     * @param x0 x0
     * @param d d
     * @param t t
     * @return [com.company.EquationSystem.discrepancy](x0 + t * d)
     */
    fun discrepancy(x0: DoubleArray, d: DoubleArray, t: Double): Double {
        val x = DoubleArray(n)
        for (i in 0 until n) {
            x[i] = x0[i] + t * d[i]
        }
        return discrepancy(x)
    }

    /** Finds solution of equation F'dx+Fx=0, where F is matrix of fi
     * @param x initial function argument
     * @return dx
     */
    fun linearDerivativeSolution(x: DoubleArray): DoubleArray {
        val b = DoubleArray(n, init = {i -> -functions[i].calculate(x)})
        val matrix = Array(n, init = {i -> functions[i].totalDerivative(x)})
        val m = Matrix(matrix)
        return m.gaussMethod(b)
    }

    /**
     * Finds local discrepancy minimum on line x(t) = x + t * direction
     * @param x initial argument
     * @param d line direction
     * @return t| x(t) = x + t * direction is local discrepancy minimum
     */
    fun localMinimum(x: DoubleArray, d: DoubleArray): Double {
        var cur = discrepancy(x)
        cur = Math.min(cur, discrepancy(x, d, 1.0))
        var r = 1.0
        var dr: Double
        do {
            r *= 2.0
            dr = discrepancy(x, d, r)
            cur = Math.min(cur, dr)
        } while (dr <= cur)
        val f = object : Function() {
            override fun calculate(arg: DoubleArray): Double {
                return discrepancy(x, d, arg[0])
            }
        }
        return gradientDescent(f, 1.0, 0.5, GRADIENT_DESCENT_PRECISION)
    }

    /**
     * @param eps precision of finding x
     * @param maxIterations maximum iterations count
     * @return argument x, discrepancy(x) < eps
     */
    fun universalMethod(eps: Double, maxIterations: Long): DoubleArray {
        val x = DoubleArray(n)
        for (i in 0 until n) {
            x[i] = random.nextDouble()
        }
        for (q in 0 until maxIterations) {
            val dx = linearDerivativeSolution(x)
            val k = localMinimum(x, dx)
            for (i in 0 until n) {
                x[i] += k * dx[i]
            }
            if (getNorm(dx) < eps) break
        }
        return x
    }

    companion object {

        private val random = RandomHolder.random
        private val GRADIENT_DESCENT_PRECISION = 1e-6


        /**
         * Finds infinity norm of vector x
         * @param x vector
         * @return ||x||{Inf}
         */
        fun getNorm(x: DoubleArray): Double {
            var norm = Math.abs(x[0])
            for (i in 1 until x.size) {
                if (Math.abs(x[i]) > norm) {
                    norm = Math.abs(x[i])
                }
            }
            return norm
        }

        /**
         * Finds argument of rough function minimum
         * @param f unary function
         * @param x0 initial argument
         * @param initialStep initial dx
         * @param precision minimum dx
         * @return x| f(x) is rough minimum
         */
        fun gradientDescent(f: Function, x0: Double, initialStep: Double, precision: Double): Double {
            val arg = DoubleArray(1)
            arg[0] = x0
            var derivative = f.totalDerivative(arg)[0]
            var step = initialStep
            var x = x0
            var x1: Double
            var min = f.calculate(arg)
            while (step > precision) {
                if (derivative < 0) {
                    x1 = x + step
                } else {
                    x1 = x - step
                }
                arg[0] = x1
                val cur = f.calculate(arg)
                if (cur < min) {
                    min = cur
                    x = x1
                    arg[0] = x
                    derivative = f.totalDerivative(arg)[0]
                } else {
                    step /= 2.0
                }
            }
            return x
        }
    }
}