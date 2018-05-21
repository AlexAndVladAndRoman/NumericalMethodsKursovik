import java.util.Random

class Matrix {

    private var n: Int = 0
    private var a: Array<DoubleArray>? = null

    private var norm = -1.0

    val matrixCopy: Array<DoubleArray>
        get() {
            val result = Array(n) { DoubleArray(n) }
            for (i in 0 until n) {
                for (j in 0 until n) {
                    result[i][j] = a!![i][j]
                }
            }
            return result
        }

    val isDiagonalDominant: Boolean
        get() {
            for (i in 0 until n) {
                var sum = 0.0
                for (j in 0 until i) {
                    sum += Math.abs(a!![i][j])
                }
                for (j in i + 1 until n) {
                    sum += Math.abs(a!![i][j])
                }
                if (sum >= Math.abs(a!![i][i])) {
                    return false
                }
            }
            return true
        }

    private class InconsistentInputException(s: String) : Exception(s)

    constructor(size: Int) {
        n = size
        a = Array(n) { DoubleArray(n) }
    }

    constructor(m: Array<DoubleArray>) {
        n = m.size
        a = m
    }

    private fun solutionsFill(b: DoubleArray, min: Int, max: Int) {
        val solutions = IntArray(n)
        for (i in solutions.indices) {
            solutions[i] = random.nextInt(max - min + 1) + min
        }
        for (i in 0 until n) {
            b[i] = 0.0
            for (j in 0 until n) {
                b[i] += a!![i][j] * solutions[j]
            }
        }
    }

    fun randomFill(b: DoubleArray, min: Int, max: Int) {
        for (i in 0 until n) {
            for (j in 0 until n) {
                a!![i][j] = (random.nextInt(max - min + 1) + min).toDouble()
            }
        }
        solutionsFill(b, min, max)
    }

    fun diagonalFill(b: DoubleArray, min: Int, max: Int) {
        for (i in 0 until n) {
            for (j in 0 until n) {
                if (i == j) {
                    a!![i][j] = (random.nextInt(max - min + 1) + min).toDouble()
                } else {
                    a!![i][j] = 0.0
                }
            }
        }
        solutionsFill(b, min, max)
    }

    fun hilbertFill(b: DoubleArray, min: Int, max: Int) {
        for (i in 0 until n) {
            for (j in 0 until n) {
                a!![i][j] = 1.0 / (i.toDouble() + j.toDouble() + 1.0)
            }
        }
        solutionsFill(b, min, max)
    }

    fun diagonalDominanceFill(b: DoubleArray, min: Int, max: Int, dominanceKoef: Int) {
        for (i in 0 until n) {
            var sum = 0.0
            for (j in 0 until n) {
                a!![i][j] = (random.nextInt(max - min + 1) + min).toDouble()
                sum += Math.abs(a!![i][j])
            }
            a!![i][i] = (1 - 2 * random.nextInt(2)) * (dominanceKoef * sum + 1)
        }
        solutionsFill(b, min, max)
    }

    fun getNorm(): Double {
        if (norm < 0) {
            for (i in 0 until n) {
                var sum = 0.0
                for (j in 0 until n) {
                    sum += Math.abs(a!![i][j])
                }
                if (norm < sum) {
                    norm = sum
                }
            }
        }
        return norm
    }

    fun transform(vector: DoubleArray): DoubleArray {
        val result = DoubleArray(n)
        for (i in 0 until n) {
            result[i] = 0.0
            for (j in 0 until n) {
                result[i] += a!![i][j] * vector[j]
            }
        }
        return result
    }

    fun transposeTransform(vector: DoubleArray): DoubleArray {
        val result = DoubleArray(n)
        for (i in 0 until n) {
            result[i] = 0.0
            for (j in 0 until n) {
                result[i] += a!![j][i] * vector[j]
            }
        }
        return result
    }

    fun g(b: DoubleArray, vector: DoubleArray): DoubleArray {
        val result = transform(vector)
        for (i in 0 until n) {
            result[i] -= b[i]
        }
        return result
    }

    @Throws(InconsistentInputException::class)
    fun jacobiMethod(b: DoubleArray, maxIterations: Long, epsilon: Double,
                     zeidelMod: Boolean, relaxation: Double, check: Boolean): DoubleArray {
        val b1 = Array(n) { DoubleArray(n) }
        val b2 = Array(n) { DoubleArray(n) }

        val d = DoubleArray(n)

        for (i in 0 until n) {
            for (j in 0 until i) {
                b2[i][j] = -a!![i][j] / a!![i][i]
            }
            if (zeidelMod) {
                for (j in i + 1 until n) {
                    b1[i][j] = -a!![i][j] / a!![i][i]
                }
            } else {
                for (j in i + 1 until n) {
                    b2[i][j] = -a!![i][j] / a!![i][i]
                }
            }
            d[i] = b[i] / a!![i][i]
        }

        val major: Double

        /** Check consistency  */
        if (!zeidelMod) {
            val q = Matrix(b2).getNorm()

            major = epsilon * (1 - q) / q
            if (check) {
                if (q >= 1) {
                    throw InconsistentInputException(String.format("Inconsistent: ||B|| = %10f >= 1\n", q))
                }
                if (!isDiagonalDominant) {
                    throw InconsistentInputException("No diagonal dominance\n")
                }
            }
        } else {
            val q1 = Matrix(b1).getNorm()
            val q2 = Matrix(b2).getNorm()
            major = if (q2 != 0.0) {
                epsilon * (1 - q1) / q2
            } else {
                epsilon
            }
            if (check) {
                if (q1 + q2 >= 1) {
                    throw InconsistentInputException("Inconsistent: ||B1|| + ||B2|| >= 1 \n")
                }
            }
        }

        var prev = 0
        val x = Array(2) { DoubleArray(n) }
        for (i in 0 until n) {
            x[prev][i] = random.nextDouble()
        }

        for (k in 0 until maxIterations) {
            val next = (prev + 1) % 2
            for (i in 0 until n) {
                var value = d[i]
                for (j in 0 until n) {
                    value += b1[i][j] * x[next][j] + b2[i][j] * x[prev][j]
                }
                x[next][i] = value
            }
            var max = Math.abs(x[next][0] - x[prev][0])
            for (i in 0 until n) {
                x[next][i] = relaxation * x[next][i] + (1 - relaxation) * x[prev][i]
                if (Math.abs(x[next][i] - x[prev][i]) > max) {
                    max = Math.abs(x[next][i] - x[prev][i])
                }
            }
            prev = next
            if (max < major) break
        }

        return x[prev]
    }

    fun gaussMethod(vector: DoubleArray): DoubleArray {
        val a = matrixCopy
        val b = vector.clone()
        val order = IntArray(n)
        for (i in 0 until n) {
            order[i] = i
        }
        for (k in 0 until n) {
            /** find pivot row and pivot column */
            var maxRow = k
            var maxColumn = k
            for (i in k until n) {
                for (j in k until n) {
                    if (Math.abs(a[i][j]) > Math.abs(a[maxRow][maxColumn])) {
                        maxRow = i
                        maxColumn = j
                    }
                }
            }

            /** swap row in A matrix  */
            val temp = a[k]
            a[k] = a[maxRow]
            a[maxRow] = temp

            /** swap column in A matrix  */
            var tmp: Double
            for (i in 0 until n) {
                tmp = a[i][maxColumn]
                a[i][maxColumn] = a[i][k]
                a[i][k] = tmp
            }

            /** swap rows in answer  */
            val itmp = order[maxColumn]
            order[maxColumn] = order[k]
            order[k] = itmp

            /** swap corresponding values in constants matrix  */
            val t = b[k]
            b[k] = b[maxRow]
            b[maxRow] = t

            /** pivot within A and B  */
            for (i in k + 1 until n) {
                val factor = a[i][k] / a[k][k]
                b[i] -= factor * b[k]
                for (j in k until n) {
                    a[i][j] -= factor * a[k][j]
                }
            }
        }

        val v = DoubleArray(n)
        for (i in n - 1 downTo 0) {
            var sum = 0.0
            for (j in i + 1 until n) {
                sum += a[i][j] * v[j]
            }
            v[i] = (b[i] - sum) / a[i][i]
        }
        val solution = DoubleArray(n)
        for (i in 0 until n) {
            solution[order[i]] = v[i]
        }
        return solution
    }

    private fun conjugateGradientsMethod(b: DoubleArray, r: Long): DoubleArray {
        var b = b
        val sym = Array(n) { DoubleArray(n) }
        for (i in 0 until n) {
            for (j in 0 until n) {
                for (k in 0 until n) {
                    sym[i][j] += a!![k][i] * a!![k][j]
                }
            }
        }
        val m = Matrix(sym)
        b = transposeTransform(b)
        val x = DoubleArray(n)
        for (i in 0 until n) {
            x[i] = random.nextDouble()
        }
        val p = m.g(b, x)
        for (i in 0 until n) {
            p[i] *= -1.0
        }
        /** With precise float operations it should give correct answer after n iterations.
         * However, we take r rounds of n iterations to be sure.  */
        for (t in 0 until n * r) {
            val Ap = m.transform(p)
            val Apxp = scalarProduct(Ap, p)
            if (Apxp == 0.0) break
            var g = m.g(b, x)
            val alpha = -scalarProduct(g, p) / Apxp
            for (i in 0 until n) {
                x[i] += alpha * p[i]
            }
            g = m.g(b, x)
            val beta = scalarProduct(Ap, g) / Apxp
            for (i in 0 until n) {
                p[i] = -g[i] + beta * p[i]
            }
        }
        return x
    }

    operator fun get(i: Int, j: Int): Double {
        return this.a!![i][j]
    }

    companion object {

        private val random = RandomHolder.random

        private fun scalarProduct(a: DoubleArray, b: DoubleArray): Double {
            var result = 0.0
            for (i in a.indices) {
                result += a[i] * b[i]
            }
            return result
        }
    }
}

