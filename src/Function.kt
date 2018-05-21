import java.util.Arrays

abstract class Function {

    abstract fun calculate(x: DoubleArray): Double

    fun totalDerivative(x: DoubleArray): DoubleArray {
        val xn = Arrays.copyOf(x, x.size)
        val y = calculate(x)
        val res = DoubleArray(x.size)
        for (i in x.indices) {
            xn[i] += EPS
            res[i] = (calculate(xn) - y) / EPS
            xn[i] = x[i]
        }
        return res
    }

    companion object {
        private val EPS = 1e-6
    }
}