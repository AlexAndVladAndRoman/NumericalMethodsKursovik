import java.io.BufferedReader
import java.io.FileReader
import java.io.IOException
import java.util.*

object DataHolder {
    val ATMOSPHERIC_PRESSURE = 100000.0
    val R = 8.314
    private val data = HashMap<String, HashMap<String, String>>()
    private var hasData = false

    @Throws(IOException::class)
    private fun readData() {
        val `in` = BufferedReader(FileReader("Bank_TD_Fragment.dat"))
        for (i in 0..5) `in`.readLine()
        var tokens = StringTokenizer(`in`.readLine())
        val parameters = ArrayList<String>()
        tokens.nextToken()
        while (tokens.hasMoreTokens()) {
            val parameter = tokens.nextToken().toUpperCase()
            parameters.add(parameter)
            data[parameter] = HashMap()
        }
        `in`.readLine()
        tokens = StringTokenizer(`in`.readLine())
        while (tokens.hasMoreTokens()) {
            val chemicalAgent = tokens.nextToken().toUpperCase()
            for (parameter in parameters) {
                data[parameter]!![chemicalAgent] = tokens.nextToken()
            }
            tokens = StringTokenizer(`in`.readLine())
        }
        val densityMap = HashMap<String, String>()
        densityMap["AL"] = "2690"
        densityMap["GA"] = "5900"
        densityMap["ALN"] = "3200"
        densityMap["GAN"] = "6150"
        data["DENSITY"] = densityMap

    }

    fun getData(parameter: String, chemicalAgent: String): String? {
        if (!hasData) {
            try {
                readData()
                hasData = true
            } catch (e: IOException) {
                e.printStackTrace()
                return null
            }

        }
        val values = data[parameter.toUpperCase()]
        return values!![chemicalAgent.toUpperCase()]
    }

    fun getDouble(parameter: String, chemicalAgent: String): Double {
        return java.lang.Double.parseDouble(getData(parameter, chemicalAgent)!!)
    }

    fun getD(chemicalAgent: String, T: Double): Double {
        val sigmaIN2 = (getDouble("sigma", chemicalAgent) + getDouble("sigma", "N2")) / 2
        val epsIN2 = Math.sqrt(getDouble("eps", chemicalAgent) * getDouble("eps", "N2"))
        val omega11 = 1.074 * Math.pow(T / epsIN2, -0.1604)
        val muIN2 = 2.0 * getDouble("mu", chemicalAgent) * getDouble("mu", "N2") / (getDouble("mu", chemicalAgent) + getDouble("mu", "N2"))
        return 0.02628 * Math.pow(T, 1.5) / (ATMOSPHERIC_PRESSURE * sigmaIN2 * omega11 * Math.sqrt(muIN2))
    }

    fun getG(chemicalAgent: String, T: Double): Double {
        val x = T / 10000
        return getDouble("H", chemicalAgent) - T * (getDouble("f1", chemicalAgent) + getDouble("f2", chemicalAgent) * Math.log(x)
                + getDouble("f3", chemicalAgent) / (x * x) + getDouble("f4", chemicalAgent) / x + getDouble("f5", chemicalAgent) * x
                + getDouble("f6", chemicalAgent) * x * x + getDouble("f7", chemicalAgent) * x * x * x)
    }

    fun getK(T: Double, number: Int): Double {
        when (number) {
            1 ->
                //2 HCl + 2 Al = 2 AlCl + H2
                return Math.exp(-(2 * getG("HCl", T) + 2 * getG("Al", T) - 2 * getG("AlCl", T) - getG("H2", T)) / (R * T)) / ATMOSPHERIC_PRESSURE
            2 ->
                // 2 HCl + Al = AlCl2 + H2
                return Math.exp(-(2 * getG("HCl", T) + getG("Al", T) - getG("AlCl2", T) - getG("H2", T)) / (R * T))
            3 ->
                // 6 HCl + 2 Al = 2 AlCl3 + 3 H2
                return Math.exp(-(6 * getG("HCl", T) + 2 * getG("Al", T) - 2 * getG("AlCl3", T) - 3 * getG("H2", T)) / (R * T)) * ATMOSPHERIC_PRESSURE
            4 ->
                // 2 HCl + 2 Ga = 2 GaCl + H2
                return Math.exp(-(2 * getG("HCl", T) + 2 * getG("Ga", T) - 2 * getG("GaCl", T) - getG("H2", T)) / (R * T)) / ATMOSPHERIC_PRESSURE
            5 ->
                // 2 HCl + Ga = GaCl2 + H2
                return Math.exp(-(2 * getG("HCl", T) + getG("Ga", T) - getG("GaCl2", T) - getG("H2", T)) / (R * T))
            6 ->
                //6 HCl + 2 Ga = 2 GaCl3 + 3 H2
                return Math.exp(-(6 * getG("HCl", T) + 2 * getG("Ga", T) - 2 * getG("GaCl3", T) - 3 * getG("H2", T)) / (R * T)) * ATMOSPHERIC_PRESSURE
            7 ->
                // AlCl + NH3 = AlN + HCl + H2
                return Math.exp(-(getG("AlCl", T) + getG("NH3", T) - getG("AlN", T) - getG("HCl", T) - getG("H2", T)) / (R * T))
            8 ->
                // 2 AlCl2 + 2 NH3 = 2 AlN + 4 HCl + H2
                return Math.exp(-(2 * getG("AlCl2", T) + 2 * getG("NH3", T) - 2 * getG("AlN", T) - 4 * getG("HCl", T) - getG("H2", T)) / (R * T)) / ATMOSPHERIC_PRESSURE
            9 ->
                // AlCl3 + NH3 = AlN + 3 HCl
                return Math.exp(-(getG("AlCl3", T) + getG("NH3", T) - getG("AlN", T) - 3 * getG("HCl", T)) / (R * T)) / ATMOSPHERIC_PRESSURE
            10 ->
                // GaCl + NH3 = GaN + HCl + H2
                return Math.exp(-(getG("GaCl", T) + getG("NH3", T) - getG("GaN", T) - getG("HCl", T) - getG("H2", T)) / (R * T))
            11 ->
                // 2 GaCl2 + 2 NH3 = 2 GaN + 4 HCl + H2
                return Math.exp(-(2 * getG("GaCl2", T) + 2 * getG("NH3", T) - 2 * getG("GaN", T) - 4 * getG("HCl", T) - getG("H2", T)) / (R * T)) / ATMOSPHERIC_PRESSURE
            12 ->
                // GaCl3 + NH3 = GaN + 3HCl
                return Math.exp(-(getG("GaCl3", T) + getG("NH3", T) - getG("GaN", T) - 3 * getG("HCl", T)) / (R * T)) / ATMOSPHERIC_PRESSURE
            else -> return 0.0
        }
    }
}