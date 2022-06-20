//
// Created by wulong on 1/30/19.
//

#include "EMAlgorithm.h"

Gaussian::Gaussian(double mu, double sigma) : m_mu(mu), m_sigma(sigma) {}

double Gaussian::getPdf(double x) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m_mu) / m_sigma;
    return inv_sqrt_2pi / m_sigma * std::exp(-0.5 * a * a);
}

void Gaussian::update(const vector<double> &sample, const vector<double> &R) {
    vector<double> dotproduct(R.size(), 0);
    for (int j = 0; j < sample.size(); j++) {
        dotproduct[j] = sample[j] * R[j];
    }
    double stable_mu = calcSum(dotproduct, false) / calcSum(R, false);

    m_mu = stable_mu;
    for (int j = 0; j < sample.size(); j++) {
        dotproduct[j] = (sample[j] - m_mu) * (sample[j] - m_mu) * R[j];
    }
    m_sigma = sqrt(calcSum(dotproduct, false) / calcSum(R, false));
}

void Gaussian::display() {
    cout << setprecision(8) << "N(" << m_mu << "," << m_sigma << ") ";
}

double Gaussian::getCdf(double x) {
    double a = (x - m_mu) / sqrt(2) / m_sigma;
    return 0.5 * (1 + erf(a));
}

Gumbel::Gumbel(double mu, double sigma) : m_mu(mu), m_sigma(sigma) {
    updateAlphaBeta();
}

void Gumbel::updateAlphaBeta() {
    const double gamma = 0.57721566490153286060651209008240243;
    const double pi = 3.141592653589793238462643383279;
    m_beta = m_sigma * sqrt(6) / pi;
    m_alpha = m_mu - m_beta * gamma;
}

double Gumbel::getPdf(double x) {
    double z = (x - m_alpha) / m_beta;
    double a = z + exp(-z);
    return 1 / m_beta * exp(-a);
}

void Gumbel::update(const vector<double> &sample, const vector<double> &R) {
    double sum = 0;
    double sumr = 0;
    // calculate mean
    for (int j = 0; j < sample.size(); j++) {
        sum += (sample[j] * R[j]);
        sumr += R[j];
    }
    m_mu = sum / sumr;

    sum = 0;
    for (int j = 0; j < sample.size(); j++) {
        sum += (sample[j] - m_mu) * (sample[j] - m_mu) * R[j];
    }
    sum /= sumr;
    m_sigma = sqrt(sum);

    updateAlphaBeta();
}

void Gumbel::display() {
    cout << "Gumbel(" << m_mu << " " << m_sigma << ") ";
}

EM::EM(vector<Distribution *> d, vector<double> &sample) : Sample(sample) {
    cout << "Model Mixture Num: " << d.size() << endl;
    D.assign(d.begin(), d.end());
    W.assign(d.size(), 1.0 / d.size());
    r.assign(d.size(), vector<double>(sample.size(), 1.0 / d.size()));
}

double EM::pdfMix(double x) {
    double sum = 0;
    for (int i = 0; i < D.size(); i++) {
        sum += W[i] * D[i]->getPdf(x);
    }
    return sum;
}

void EM::display() {
    cout << "Model: ";
    for (int i = 0; i < W.size(); i++) {
        if (i > 0) cout << " + ";
        cout << W[i] << " x ";
        D[i]->display();
    }
    cout << endl;
}

void EM::updateR() {
    for (int i = 0; i < D.size(); i++) {
        for (int j = 0; j < Sample.size(); j++) {
            r[i][j] = W[i] * D[i]->getPdf(Sample[j]);
        }
    }

    // normalization
    for (int j = 0; j < Sample.size(); j++) {
        double s = 0;
        for (int i = 0; i < D.size(); i++) {
            s += r[i][j];
        }

        for (int i = 0; i < D.size(); i++) {
            r[i][j] /= s;
        }
    }
}

void EM::updateModel() {
    for (int i = 0; i < D.size(); i++) {
        double stable_s = calcSum(r[i], false);
        W[i] = stable_s / Sample.size();
    }

    // update model
    for (int i = 0; i < D.size(); i++) {
        D[i]->update(Sample, r[i]);
    }
}

double EM::loglikelihood() // todo: should be less than ZERO to be fixed
{
    double s = 0;
    for (int j = 0; j < Sample.size(); j++) {
        s += log(pdfMix(Sample[j]));
    }
    return s;
}

void EM::run(int n) {
    cout << "EM likelihood: " << endl;
    while (n--) {
        cout << " - - " << loglikelihood() << flush;
        updateR();
        updateModel();
    }
    double likelihood = loglikelihood();
    cout << "EM likelihood " << likelihood << endl;
    display();
}

double calcSum(const vector<double> &x, bool stable) {
    double s = 0;
    if (not stable) {
        s = accumulate(x.begin(), x.end(), 0.0L);
    } else {
        vector<double> y(x.begin(), x.end());
        std::sort(y.begin(), y.end());

        for (int i = 0; i < y.size(); i++) {
            s += y[i];
        }
    }
    return s;
}

void visualizeHistogram(const vector<double> &s) {
    int N = s.size();
    double maxval = *std::max_element(s.begin(), s.end());
    double minval = *std::min_element(s.begin(), s.end());
    // plot figure
    int totalDots = 60;
    int countsPerDot = N / totalDots;
    int binNum = 20;

    std::map<double, int> hist{};
    for (int i = 0; i < N; ++i) {
        ++hist[std::round((s[i] - minval) / (maxval - minval) * binNum)];
    }
    for (auto p : hist) {
        std::cout << std::setw(9) << setprecision(6)
                  << p.first * (maxval - minval) / binNum + minval << ' '
                  << std::string(p.second / (countsPerDot), '*')
                  << "|" << p.second << '\n';
    }
}
