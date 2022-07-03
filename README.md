// This source code is subject to the terms of the Mozilla Public License 2.0 at https://mozilla.org/MPL/2.0/
// © loxx

//@version=5
indicator("Adaptive, Jurik-Filtered, Floating RSI [Loxx]", shorttitle="AJFFRSI [Loxx]", timeframe="", overlay = false, timeframe_gaps=true, max_bars_back = 3000)

greencolor = #2DD204  
redcolor = #D2042D 

RMA(x, t) =>
    EMA1 = x
    EMA1 := na(EMA1[1]) ? x : (x - nz(EMA1[1])) * (1/t) + nz(EMA1[1])
    EMA1

_f_hp(_src, max_len) =>
    var c = 360 * math.pi / 180
    _alpha = (1 - math.sin(c / max_len)) / math.cos(c / max_len)
    _hp = 0.0
    _hp := 0.5 * (1 + _alpha) * (_src - nz(_src[1])) + _alpha * nz(_hp[1])
    _hp

_f_ess(_src, _len) =>
    var s = 1.414
    _a = math.exp(-s * math.pi / _len)
    _b = 2 * _a * math.cos(s *  math.pi  / _len)
    _c2 = _b
    _c3 = -_a * _a
    _c1 = 1 - _c2 - _c3
    _out = 0.0
    _out := _c1 * (_src + nz(_src[1])) / 2 + _c2 * nz(_out[1], nz(_src[1], _src)) + _c3 * nz(_out[2], nz(_src[2], nz(_src[1], _src)))
    _out

_auto_dom(src, min_len, max_len, ave_len) =>
    var c = 2 * math.pi
    var s = 1.414
    avglen = ave_len

    filt = _f_ess(_f_hp(src, max_len), min_len)
    arr_size = max_len * 2
    
    var corr = array.new_float(arr_size, initial_value=0)
    var cospart = array.new_float(arr_size, initial_value=0)
    var sinpart = array.new_float(arr_size, initial_value=0)
    var sqsum = array.new_float(arr_size, initial_value=0)
    var r1 = array.new_float(arr_size, initial_value=0)
    var r2 = array.new_float(arr_size, initial_value=0)
    var pwr = array.new_float(arr_size, initial_value=0)
    
    for lag = 0 to max_len by 1
        m = avglen == 0 ? lag : avglen
        Sx = 0.0, Sy = 0.0,  Sxx = 0.0, Syy = 0.0, Sxy = 0.0
        for i = 0 to m - 1 by 1
            x = nz(filt[i])
            y = nz(filt[lag + i])
            Sx += x
            Sy += y
            Sxx += x * x
            Sxy += x * y
            Syy += y * y
            Syy
        if (m * Sxx - Sx * Sx) * (m * Syy - Sy * Sy) > 0
            array.set(corr, lag, (m * Sxy - Sx * Sy) / math.sqrt((m * Sxx - Sx * Sx) * (m * Syy - Sy * Sy)))
    
    for period = min_len to max_len by 1
        array.set(cospart, period, 0)
        array.set(sinpart, period, 0)
        for n = ave_len to max_len by 1
            array.set(cospart, period, nz(array.get(cospart, period)) + nz(array.get(corr, n)) * math.cos(c * n / period))
            array.set(sinpart, period, nz(array.get(sinpart, period)) + nz(array.get(corr, n)) * math.sin(c * n / period))
        array.set(sqsum, period, math.pow(nz(array.get(cospart, period)), 2) + math.pow(nz(array.get(sinpart, period)), 2))
    
    for period = min_len to max_len by 1
        array.set(r2, period, nz(array.get(r1, period)))
        array.set(r1, period, 0.2 * math.pow(nz(array.get(sqsum, period)), 2) + 0.8 * nz(array.get(r2, period)))

    maxpwr = 0.0
    
    for period = min_len to max_len by 1
        if nz(array.get(r1, period)) > maxpwr
            maxpwr := nz(array.get(r1, period))

    for period = ave_len to max_len by 1
        array.set(pwr, period, nz(array.get(r1, period)) / maxpwr)
    
    dominantcycle = 0.0, peakpwr = 0.0
    
    for period = min_len to max_len by 1
        if nz(array.get(pwr, period)) > peakpwr 
            peakpwr := nz(array.get(pwr, period))

    spx = 0.0, sp = 0.0

    for period = min_len to max_len by 1
        if peakpwr >= 0.25 and nz(array.get(pwr, period)) >= 0.25
            spx += period * nz(array.get(pwr, period))
            sp += nz(array.get(pwr, period))

    dominantcycle := sp != 0 ? spx / sp : dominantcycle
    dominantcycle := sp < 0.25  ?  dominantcycle[1] : dominantcycle
    dominantcycle := dominantcycle < 1 ? 1 : dominantcycle 
    dom_in = math.min(math.max(dominantcycle, min_len), max_len)
    dom_in

jma(src, len, phase, power) =>
    phaseRatio = phase < -100 ? 0.5 : phase > 100 ? 2.5 : phase / 100 + 1.5
    beta = 0.45 * (len - 1) / (0.45 * (len - 1) + 2)
    alpha = math.pow(beta, power)
    jma1 = 0.0
    e0 = 0.0
    e0 := (1 - alpha) * src + alpha * nz(e0[1])
    e1 = 0.0
    e1 := (src - e0) * (1 - beta) + beta * nz(e1[1])
    e2 = 0.0
    e2 := (e0 + phaseRatio * e1 - nz(jma1[1])) * math.pow(1 - alpha, 2) + math.pow(alpha, 2) * nz(e2[1])
    jma1 := e2 + nz(jma1[1])
    jma1

volty(src, len)=>
    len1 = math.max(math.log(math.sqrt(0.5 * (len-1))) / math.log(2.0) + 2.0,0)
    pow1 = math.max(len1 - 2.0, 0.5)
    
    voltya = 2.0, avolty = 4.0, vsum = 3.0, bsmax = 0.0, bsmin = 1.0, avgLen = 65
    
    del1 = src - nz(bsmax[1])
    del2 = src -  nz(bsmin[1])
    
    if(math.abs(del1) > math.abs(del2)) 
        voltya := math.abs(del1) 
        
    if(math.abs(del1) < math.abs(del2))
        voltya := math.abs(del2) 
        
    vsum := nz(vsum[1]) + 0.1 * (voltya - nz(voltya[9]))
    avg = ta.sma(vsum, avgLen)
    
    avolty := avg                                           
    dVolty = avolty > 0  ? voltya / avolty : 0
    
    if (dVolty > math.pow(len1, 1.0/pow1))
        dVolty := math.pow(len1, 1.0/pow1)
        
    if (dVolty < 1)                      
        dVolty := 1.0
    
    pow2 = math.pow(dVolty, pow1)
    len2 = math.sqrt(0.5 * (len - 1)) * len1
    Kv = math.pow(len2 / (len2 + 1), math.sqrt(pow2))
    
    bsmax := del1 > 0 ? src : src - Kv * del1
    bsmin := del2 < 0 ? src : src - Kv * del2
    pow2
_rsx_rsi(src, len)=>
    f8 = 100 * src
    f10 = nz(f8[1])
    v8 = f8 - f10
    f18 = 3 / (len + 2)
    f20 = 1 - f18
    f28 = 0.0
    f28 := f20 * nz(f28[1]) + f18 * v8
    f30 = 0.0
    f30 := f18 * f28 + f20 * nz(f30[1])
    vC = f28 * 1.5 - f30 * 0.5
    f38 = 0.0
    f38 := f20 * nz(f38[1]) + f18 * vC
    f40 = 0.0
    f40 := f18 * f38 + f20 * nz(f40[1])
    v10 = f38 * 1.5 - f40 * 0.5
    f48 = 0.0
    f48 := f20 * nz(f48[1]) + f18 * v10
    f50 = 0.0
    f50 := f18 * f48 + f20 * nz(f50[1])
    v14 = f48 * 1.5 - f50 * 0.5
    f58 = 0.0
    f58 := f20 * nz(f58[1]) + f18 * math.abs(v8)
    f60 = 0.0
    f60 := f18 * f58 + f20 * nz(f60[1])
    v18 = f58 * 1.5 - f60 * 0.5
    f68 = 0.0
    f68 := f20 * nz(f68[1]) + f18 * v18
    f70 = 0.0
    f70 := f18 * f68 + f20 * nz(f70[1])
    v1C = f68 * 1.5 - f70 * 0.5
    f78 = 0.0
    f78 := f20 * nz(f78[1]) + f18 * v1C
    f80 = 0.0
    f80 := f18 * f78 + f20 * nz(f80[1])
    v20 = f78 * 1.5 - f80 * 0.5
    f88_ = 0.0
    f90_ = 0.0
    f88 = 0.0
    f90_ := nz(f90_[1]) == 0 ? 1 : nz(f88[1]) <= nz(f90_[1]) ? nz(f88[1]) + 1 : nz(f90_[1]) + 1
    f88 := nz(f90_[1]) == 0 and len - 1 >= 5 ? len - 1 : 5
    f0 = f88 >= f90_ and f8 != f10 ? 1 : 0
    f90 = f88 == f90_ and f0 == 0 ? 0 : f90_
    v4_ = f88 < f90 and v20 > 0 ? (v14 / v20 + 1) * 50 : 50
    rsx = v4_ > 100 ? 100 : v4_ < 0 ? 0 : v4_
    rsx

_wilders_rsi(x, y) => 
    u = math.max(x - x[1], 0) 
    d = math.max(x[1] - x, 0)  
    rs = RMA(u, y) / RMA(d, y)
    res = 100 - 100 / (1 + rs)
    res
    
_rapid_rsi(src, len)=>
    upSum = math.sum(math.max(ta.change(src), 0), len)
    dnSum = math.sum(math.max(-ta.change(src), 0), len)
    rrsi = dnSum == 0 ? 100 : upSum == 0 ? 0 : 100 - 100 / (1 + upSum / dnSum)
    rrsi
    
adapt = input.string("Fixed", "Calculation type", options = ["Fixed", "Autocorrelation Adaptive"], group = "Basic Settings")

rsiType = input.string("Wilders", "RSI Type", options =["Wilders", "RSX", "Rapid"], group = "Basic Settings")
src = input.source(close, "Source", group = "Basic Settings")
len = input.int(15, "Length", group = "Basic Settings")
MinMaxPeriod = input.int(100, "Max floating period", group = "Basic Settings")
phs = input.int(50, title= "Jurik Phase", group = "Basic Settings")

overb = input.int(80, "Overbought", group = "Thresholds")
overs = input.int(20, "Oversold", group = "Thresholds")

auto_src = input.source(close, title='APA Source', group = "Autocorrelation Periodgram Algorithm")
auto_min = input.int(8, minval = 1, title='APA Min Length', group = "Autocorrelation Periodgram Algorithm")
auto_max = input.int(48, minval = 1, title='APA Max Length', group = "Autocorrelation Periodgram Algorithm")
auto_avg = input.int(3, minval = 1, title='APA Average Length', group = "Autocorrelation Periodgram Algorithm")

colorbars = input.bool(false, "Color bars?", group = "UI Options")

auto = _auto_dom(auto_src, auto_min, auto_max, auto_avg)
float out = adapt == "Fixed" ? len : int(auto) < 1 ? 1 : int(auto) 

rsi = rsiType == "Wilders" ? _wilders_rsi(src, int(auto)) : rsiType == "Rapid" ? _rapid_rsi(src, int(auto)) :  _rsx_rsi(src, int(auto))

volty = volty(rsi, out)
rsi := jma(rsi, out, phs, volty)

up = overb 
dn = overs

min = ta.lowest(rsi, 100)
max = ta.highest(rsi, 100)

rng = max - min

levelDn = min + rng * dn / 100.0
levelMi = min + rng * 0.5
levelUp = min + rng * up / 100.0

rsiplot = plot(rsi, color = color.white, linewidth = 2)
dnplot = plot(levelDn, color = greencolor, style = plot.style_line)
midplot = plot(levelMi, color = color.gray, style = plot.style_line)
upplot = plot(levelUp, color = redcolor, style = plot.style_line)


fill(rsiplot, upplot, color = rsi >= levelUp ? redcolor :na )
fill(dnplot, rsiplot, color = rsi <= levelDn ? greencolor :na )
barcolor(colorbars ? rsi > levelMi ? greencolor : redcolor : na)










    
