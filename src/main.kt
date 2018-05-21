import java.io.*
import java.util.*
import com.sun.xml.internal.ws.streaming.XMLStreamReaderUtil.close
import java.io.PrintWriter
import java.util.Locale
import java.util.HashMap
import java.io.IOException


// We allow answer pressures in [-ALLOWED_DISCREPANCY; ATMOSPHERIC_PRESSURE + ALLOWED_DISCREPANCY]
val ALLOWED_DISCREPANCY = 1000.0

fun solveAlClx(pressure: Map<String, Double>, T: Double, delta: Double, out: PrintWriter) {
    val chemicalAgent = arrayOf("HCl", "AlCl", "AlCl2", "AlCl3", "H2")
    val k1 = DataHolder.getK(T, 1)
    val k2 = DataHolder.getK(T, 2)
    val k3 = DataHolder.getK(T, 3)
    // x[i] = Pe(chemicalAgent[i])
    val functions = arrayOf(
            // 2 HCl + 2 Al = 2 AlCl + H2
            // Pe(HCl)^2 = K1 * Pe(AlCl)^2 * Pe(H2)
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return x[0] * x[0] - k1 * x[1] * x[1] * x[4]
            }
        },
            // 2 HCl + Al = AlCl2 + H2
            // Pe(HCl) ^ 2 = K2 * Pe(AlCl2) * Pe(H2)
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return x[0] * x[0] - k2 * x[2] * x[4]
            }
        },
            // 6 HCl + 2 Al = 2 AlCl3 + 3 H2
            // Pe(HCl)^6 = K3 * Pe(AlCl3)^2 * Pe(H2)^3
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return Math.pow(x[0], 6.0) - k3 * x[3] * x[3] * x[4] * x[4] * x[4]
            }
        })

    val p = DoubleArray(5)
    val d = DoubleArray(5)
    for (i in 0..4) {
        p[i] = pressure.get(chemicalAgent[i])!!
        d[i] = DataHolder.getD(chemicalAgent[i], T)
    }
    // G(H) = G(HCl) + 2 * G(H2) = 0
    // D(HCl) * (Pg(HCl) - Pe(HCl)) + 2 * D(H2) * (Pg(H2) - Pe(H2)) = 0
    functions[3] = object : Function() {
        override fun calculate(x: DoubleArray): Double {
            return d[0] * (p[0] - x[0]) + 2.0 * d[4] * (p[4] - x[4])
        }
    }
    // G(Cl) = G(HCl) + G(AlCl) + 2 * G(AlCl2) + 3 * G(AlCl3) = 0
    // D(HCl) * (Pg(HCl) - Pe(HCl)) + D(AlCl) * (Pg(AlCl) - Pe(AlCl)) + 2 * D(AlCl2) * (Pg(AlCl2) - Pe(AlCl2))
    // + 3 * D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) = 0
    functions[4] = object : Function() {
        override fun calculate(x: DoubleArray): Double {
            return d[1] * (p[1] - x[1]) + 2.0 * d[2] * (p[2] - x[2]) + 3.0 * d[3] * (p[3] - x[3]) + d[0] * (p[0] - x[0])
        }
    }
    val equationSystem = EquationSystem(functions)
    var correct = false
    var x: DoubleArray? = null
    while (!correct) {
        x = equationSystem.universalMethod(1e-12, 1000000)
        correct = true
        for (i in 0..4) {
            correct = correct and (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.ATMOSPHERIC_PRESSURE + ALLOWED_DISCREPANCY)
        }
    }
    println("T = $T")
    out.println("T = $T")
    for (i in 0..4) {
        println("Pe(" + chemicalAgent[i] + ") = " + x!![i])
        out.println("Pe(" + chemicalAgent[i] + ") = " + x[i])
    }
    val g = DoubleArray(5)
    for (i in 0..4) {
        g[i] = d[i] * (p[i] - x!![i]) / (8314.0 * T * delta)
        println("G(" + chemicalAgent[i] + ") = " + g[i])
        out.println("G(" + chemicalAgent[i] + ") = " + g[i])
    }
    val v = (g[1] + g[2] + g[3]) * (DataHolder.getDouble("mu", "Al") / DataHolder.getDouble("density", "Al")) * 1000000000.0
    println("Ve(Al) = $v")
    out.println("Ve(Al) = $v")
}

fun solveGaClx(pressure: Map<String, Double>, T: Double, delta: Double, out: PrintWriter) {
    val chemicalAgent = arrayOf("HCl", "GaCl", "GaCl2", "GaCl3", "H2")
    val k4 = DataHolder.getK(T, 4)
    val k5 = DataHolder.getK(T, 5)
    val k6 = DataHolder.getK(T, 6)
    // x[i] = Pe(chemicalAgent[i])
    val functions = arrayOf(
            // 2 HCl + 2 Ga = 2 GaCl + H2
            // Pe(HCl)^2 = K4 * Pe(GaCl)^2 * Pe(H2)
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return x[0] * x[0] - k4 * x[1] * x[1] * x[4]
            }
        },
            // 2 HCl + Ga = GaCl2 + H2
            // Pe(HCl)^2 = K5 * Pe(GaCl2) * Pe(H2)
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return x[0] * x[0] - k5 * x[2] * x[4]
            }
        },
            // 6 HCl + 2 Ga = 2 GaCl3 + 3 H2
            // Pe(HCl)^6 = K6 * Pe(GaCl3)^2 * Pe(H2)^3
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return Math.pow(x[0], 6.0) - k6 * x[3] * x[3] * x[4] * x[4] * x[4]
            }
        }
    )

    val p = DoubleArray(5)
    val d = DoubleArray(5)
    for (i in 0..4) {
        p[i] = pressure.get(chemicalAgent[i])!!
        d[i] = DataHolder.getD(chemicalAgent[i], T)
    }
    // G(H) = G(HCl) + 2 * G(H2) = 0
    // D(HCl) * (Pg(HCl) - Pe(HCl)) + 2 * D(H2) * (Pg(H2) - Pe(H2)) = 0
    functions[3] = object : Function() {
        override fun calculate(x: DoubleArray): Double {
            return d[0] * (p[0] - x[0]) + 2.0 * d[4] * (p[4] - x[4])
        }
    }
    // G(Cl) = G(HCl) + G(GaCl) + 2 * G(GaCl2) + 3 * G(GaCl3) = 0
    // D(HCl) * (Pg(HCl) - Pe(HCl)) + D(GaCl) * (Pg(GaCl) - Pe(GaCl)) + 2 * D(GaCl2) * (Pg(GaCl2) - Pe(GaCl2))
    // + 3 * D(GaCl3) * (Pg(GaCl3) - Pe(GaCl3)) = 0
    functions[4] = object : Function() {
        override fun calculate(x: DoubleArray): Double {
            return d[1] * (p[1] - x[1]) + 2.0 * d[2] * (p[2] - x[2]) + 3.0 * d[3] * (p[3] - x[3]) + d[0] * (p[0] - x[0])
        }
    }
    val equationSystem = EquationSystem(functions)
    var correct = false
    var x: DoubleArray? = null
    while (!correct) {
        x = equationSystem.universalMethod(1e-12, 1000000)
        correct = true
        for (i in 0..4) {
            correct = correct and (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.ATMOSPHERIC_PRESSURE + ALLOWED_DISCREPANCY)
        }
    }
    println("T = $T")
    out.println("T = $T")
    for (i in 0..4) {
        println("Pe(" + chemicalAgent[i] + ") = " + x!![i])
        out.println("Pe(" + chemicalAgent[i] + ") = " + x[i])
    }
    val g = DoubleArray(5)
    for (i in 0..4) {
        g[i] = d[i] * (p[i] - x!![i]) / (8314.0 * T * delta)
        println("G(" + chemicalAgent[i] + ") = " + g[i])
        out.println("G(" + chemicalAgent[i] + ") = " + g[i])
    }
    val v = (g[1] + g[2] + g[3]) * (DataHolder.getDouble("mu", "Ga") / DataHolder.getDouble("density", "Ga")) * 1000000000.0
    println("Ve(Ga) = $v")
    out.println("Ve(Ga) = $v")
}

fun solveAlGaN(pressure: Map<String, Double>, T: Double, delta: Double, out: PrintWriter) {
    val chemicalAgent = arrayOf("HCl", "GaCl", "NH3", "AlCl3", "H2")
    val p = DoubleArray(5)
    val d = DoubleArray(5)
    for (i in 0..4) {
        p[i] = pressure.get(chemicalAgent[i])!!
        d[i] = DataHolder.getD(chemicalAgent[i], T)
    }
    val k9 = DataHolder.getK(T, 9)
    val k10 = DataHolder.getK(T, 10)
    // x[i] = Pe(chemicalAgent[i]), i = 0..4
    // x[5] = x = G(AlCl3) / (G(AlCl3) + G(GaCl))
    val functions = arrayOf(
        // AlCl3 + NH3 = AlN + 3 HCl
        // Pe(AlCl3) * Pe(NH3) = K9 * x * Pe(HCl)^3
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return x[3] * x[2] - k9 * x[5] * x[0] * x[0] * x[0]
            }
        },
        // GaCl + NH3 = GaN + HCl + H2
        // Pe(GaCl) * Pe(NH3) = K10 * (1 - x) * Pe(HCl) * Pe(H2)
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return d[0] * (p[0] - x[0]) + 2.0 * d[4] * (p[4] - x[4]) + 3.0 * d[2] * (p[2] - x[2])
            }
        },
        // G(H) = G(HCl) + 2 * G(H2) + 3 * G(NH3) = 0
        // D(HCl) * (Pg(HCl) - Pe(HCl)) + 2 * D(H2) * (Pg(H2) - Pe(H2)) + 3 * D(NH3) * (Pg(NH3) - Pe(NH3))
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return x[1] * x[2] - k10 * (1 - x[5]) * x[0] * x[4]
            }
        },
        // G(Cl) = 3 * G(AlCl3) + G(GaCl) + G(HCl) = 0
        // 3 * D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) + D(GaCl) * (Pg(GaCl) - Pe(GaCl)) + D(HCl) * (Pg(HCl) - Pe(HCl)) = 0
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return 3.0 * d[3] * (p[3] - x[3]) + d[1] * (p[1] - x[1]) + d[0] * (p[0] - x[0])
            }
        },
        // G(Al) + G(Ga) = G(AlCl3) + G(GaCl) = G(NH3) = G(N)
        // D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) + D(GaCl) * (Pg(GaCl) - Pe(GaCl)) = D(NH3) * (Pg(NH3) - Pe(NH3))
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return d[3] * (p[3] - x[3]) + d[1] * (p[1] - x[1]) - d[2] * (p[2] - x[2])
            }
        },
        // G(AlCl3) = x * (G(AlCl3) + G(GaCl))
        // D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) = x * (D(AlCl3) * (Pg(AlCl3) - Pe(AlCl3)) + D(GaCl) * (Pg(GaCl) - Pe(GaCl)))
        object : Function() {
            override fun calculate(x: DoubleArray): Double {
                return d[3] * (p[3] - x[3]) - x[5] * (d[1] * (p[1] - x[1]) + d[3] * (p[3] - x[3]))
            }
        }
    )
    val equationSystem = EquationSystem(functions)
    var correct = false
    var x: DoubleArray? = null
    while (!correct) {
        x = equationSystem.universalMethod(1e-12, 1000000)
        correct = true
        for (i in 0..4) {
            correct = correct and (x[i] >= -ALLOWED_DISCREPANCY && x[i] <= DataHolder.ATMOSPHERIC_PRESSURE + ALLOWED_DISCREPANCY)
        }
        correct = correct and (x[5] >= 0 && x[5] <= 1)
    }
    println("Pg(AlCl3) = " + pressure["AlCl3"])
    out.println("Pg(AlCl3) = " + pressure["AlCl3"])
    for (i in 0..4) {
        println("Pe(" + chemicalAgent[i] + ") = " + x!![i])
        out.println("Pe(" + chemicalAgent[i] + ") = " + x[i])
    }
    println("x = " + x!![5])
    out.println("x = " + x[5])
    val g = DoubleArray(5)
    for (i in 0..4) {
        g[i] = d[i] * (p[i] - x[i]) / (8314.0 * T * delta)
        println("G(" + chemicalAgent[i] + ") = " + g[i])
        out.println("G(" + chemicalAgent[i] + ") = " + g[i])
    }
    val v = (g[3] * (DataHolder.getDouble("mu", "AlN") / DataHolder.getDouble("density", "AlN")) + g[1] * (DataHolder.getDouble("mu", "GaN") / DataHolder.getDouble("density", "GaN"))) * 1000000000
    println("Vg(AlGaN) = $v")
    out.println("Vg(AlGaN) = $v")
}

@Throws(IOException::class)
fun parseFile(fileName: String): Map<String, List<Double>> {
    val result = HashMap<String, MutableList<Double>>()
    val buff = BufferedReader(FileReader(fileName))
   // var s: String
    for (line in buff.lines()) {
        val strings = line.split(" ".toRegex()).dropLastWhile { it.isEmpty() }.toTypedArray()
        if (result[strings[0]] == null) {
            result[strings[0]] = ArrayList()
        }
        result[strings[0]]!!.add(java.lang.Double.parseDouble(strings[2]))
    }
    buff.close()
    return result
}

@Throws(IOException::class)
fun main(args: Array<String>) {
    var out: PrintWriter
    val pressure = HashMap<String, Double>()
    pressure["HCl"] = 10000.0
    pressure["N2"] = 90000.0
    pressure["AlCl"] = 0.0
    pressure["AlCl2"] = 0.0
    pressure["AlCl3"] = 0.0
    pressure["GaCl"] = 0.0
    pressure["GaCl2"] = 0.0
    pressure["GaCl3"] = 0.0
    pressure["H2"] = 0.0

    // Task 1

    out = PrintWriter("task1.out")
    for (i in 35..65) {
        val T = (10 * i + 273).toDouble()
        solveAlClx(pressure, T, 0.01, out)
    }
    out.close()

    // Task 2

    out = PrintWriter("task2.out")
    for (i in 65..95) {
        val T = (10 * i + 273).toDouble()
        solveGaClx(pressure, T, 0.01, out)
    }
    out.close()

    // Task 3

    out = PrintWriter("task3_pure_N2.out")
    pressure["NH3"] = 1500.0
    pressure["HCl"] = 0.0
    println("Pure N2")
    pressure["N2"] = 98470.0
    pressure["H2"] = 0.0
    for (i in 0..30) {
        pressure["AlCl3"] = i.toDouble()
        pressure["GaCl"] = (30 - i).toDouble()
        solveAlGaN(pressure, (1100 + 273).toDouble(), 0.01, out)
    }
    out.close()
    out = PrintWriter("task3_N2_H2.out")
    println("H2/N2 = 1/9")
    pressure["N2"] = 88623.0
    pressure["H2"] = 9847.0
    for (i in 0..30) {
        pressure["AlCl3"] = i.toDouble()
        pressure["GaCl"] = (30 - i).toDouble()
        solveAlGaN(pressure, (1100 + 273).toDouble(), 0.01, out)
    }
    out.close()

    // Parsing results
    val locale = Locale("ru")
    var resMap = parseFile("task1.out")
    var t = resMap["T"]
    var valueNames = arrayOf("G(AlCl)", "G(AlCl2)", "G(AlCl3)", "Ve(Al)")
    for (valueName in valueNames) {
        out = PrintWriter("task1_$valueName.out")
        val values = resMap[valueName]
        for (i in t!!.indices) {
            out.printf(locale, "%f\t%f\n", 1 / t.get(i), Math.log(Math.abs(values!!.get(i))))
        }
        out.close()
    }
    resMap = parseFile("task2.out")
    t = resMap["T"]
    valueNames = arrayOf("G(GaCl)", "G(GaCl2)", "G(GaCl3)", "Ve(Ga)")
    for (valueName in valueNames) {
        out = PrintWriter("task2_$valueName.out")
        val values = resMap[valueName]
        for (i in t!!.indices) {
            out.printf(locale, "%f\t%f\n", 1 / t.get(i), Math.log(Math.abs(values!!.get(i))))
        }
        out.close()
    }
    resMap = parseFile("task3_pure_N2.out")
    var PgAlCl3 = resMap["Pg(AlCl3)"]
    valueNames = arrayOf("G(AlCl3)", "G(GaCl)", "Vg(AlGaN)", "x")
    for (valueName in valueNames) {
        out = PrintWriter("task3_" + valueName + "N2.out")
        val values = resMap[valueName]
        for (i in PgAlCl3!!.indices) {
            out.printf(locale, "%f\t%e\n", PgAlCl3.get(i) / 30, values!!.get(i))
        }
        out.close()
    }
    resMap = parseFile("task3_N2_H2.out")
    PgAlCl3 = resMap["Pg(AlCl3)"]
    valueNames = arrayOf("G(AlCl3)", "G(GaCl)", "Vg(AlGaN)", "x")
    for (valueName in valueNames) {
        out = PrintWriter("task3_" + valueName + "N2_H2.out")
        val values = resMap[valueName]
        for (i in PgAlCl3!!.indices) {
            out.printf(locale, "%f\t%e\n", PgAlCl3.get(i) / 30, values!!.get(i))
        }
        out.close()
    }
}
