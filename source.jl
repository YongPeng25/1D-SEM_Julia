
function source_time(t,hdur)
    decay_rate = 2.628
    PI = 3.141592653589793
    alpha = decay_rate/hdur

#    Gaussian
    # wavelet = alpha*exp(-alpha*alpha*t*t)/sqrt(PI)
#    Ricker wavelet
    wavelet = -2.0*(alpha^3)*t*exp(-alpha*alpha*t*t)/sqrt(PI)
    # wavelet = sin(t)
    wavelet
end
