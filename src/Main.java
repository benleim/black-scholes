public class Main {

    public static void main(String[] args) {

    }

    // Black-Scholes Vanilla Call
    // s = underlying price, k = strike price, r = risk-free rate, t = time (years), v = volatility
    public static double callPrice(double s, double k, double r, double t, double v) {
        double d1 = (Math.log(s/k) + (r + v*v/2) * t) / (v * Math.sqrt(t));
        double d2 = d1 - v * Math.sqrt(t);
        return s * normalCDF(d1) - (k * Math.pow(Math.E, -r*t) * normalCDF(d2));
    }

    // Black-Scholes Vanilla Put
    // s = underlying price, k = strike price, r = risk-free rate, t = time (years), v = volatility
    public static double putPrice(double s, double k, double r, double t, double v) {
        double d1 = (Math.log(s/k) + (r + v*v/2) * t) / (v * Math.sqrt(t));
        double d2 = d1 - v * Math.sqrt(t);
        double kert = k * Math.pow(Math.E, -r*t);
        return kert - s + (s*normalCDF(d1) - kert*normalCDF(d2));
    }

    // Approximation for cumulative distribution function for the standard normal distribution
    // [www.quantstart.com]
    public static double normalCDF(double x) {
        double k = 1.0/(1.0 + 0.2316419*x);
        double kSum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

        if (x >= 0.0) {
            return (1.0 - (1.0/(Math.pow(2*Math.PI,0.5)))*Math.exp(-0.5*x*x) * kSum);
        } else {
            return 1.0 - normalCDF(-x);
        }
    }
}
