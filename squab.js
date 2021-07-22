/*
    squab => small PIDgen v0.1 by @nxnine
    derived from PIDgen.c [ https://www.staff.uni-mainz.de/pommeren/PID/ ]
    derived from "An Optimal Code For Patient Identifiers" by Andreas Faldum & Klaus Pommerening
    modified for a specific use case (e.g rndwidth=0)
    XOR over mult30 to reduce collisions
*/

var sigma = "0123456789ACDEFGHJKLMNPQRTUVWXYZ";
var NN = 1073741824;
var NN_1 = 1073741823;

// init rnd_ vars
var rnd_width = 0;
var rnd_fact = 1 << (30 - rnd_width);
var rnd_lim = 1 << rnd_width;

// init keys
var k1 = 1000000000 | 1;
var k2 = 567345123 | 1;
var k3 = 765321677 | 1;
var k4 = ((k1 ^ k2) ^ k3) & NN_1;
//

// randomize a number
var rndext = function(x){
    x = x & (rnd_fact -1);
    var r = ((Math.floor(rnd_lim * Math.random()))*rnd_fact) + x;
    return r & NN_1;
};

// Rotate x cyclically by 6 bits to the right
var rot30_6 = function(x){ return ((x>>6)&16777215) | ((x&63)<<24); };
// nonlinear transform of x
var NLmix = function(x){ return (x & 1073741760) | (((x&63) + (((x>>6)&63)*((x>>24)&63)) + (((x>>12)&63)*((x>>18)&63))) & 63); };
// xor+bitmask = no collisions at 1,000,000
var xor30 = function(x, y){ return (x^y) & NN_1; };
// "encrypt" x
var encr = function(x){
    x = x & NN_1;
    // r1
    var y = rot30_6(x);
    var w = NLmix(y);
    var z = xor30(k1,w);
    // r2
    y = rot30_6(z);
    w = NLmix(y);
    z = xor30(k2,w);
    // r3
    y = rot30_6(z);
    w = NLmix(y);
    z = xor30(k3,w);
    // r4
    y = rot30_6(z);
    w = NLmix(y);
    z = xor30(k4,w);
    //
    return z;
};

// galois field arithmetic
var mult0f32 = function(x,e){
    x = x & 31;
    if (e == 0) return x;
    if (e > 3)  return x;
    var u = x >> (5-e);
    var v = u << 2;
    var w = (x << e) & 31;
    var s = (u ^ v) ^ w;
    return s;
};
var multf32 = function(x,e){
    x = x & 31;
    while (e >= 4) {
        x = mult0f32(mult0f32(x,2),2);
        e = e-4;
    }
    x = mult0f32(x,e);
    return x;
};
// Output the weighted sum
var wsum1 = function(p){
    var s = 0;
    for (var i=0; i<=5; i++) s = s ^ multf32(p[i],i+1);
    return s & 31;
};
var wsum2 = function(p){
    var s = 0;
    for (var i=0; i<=5; i++) s = s ^ multf32(p[i],2*i+2);
    return s & 31;
};
// Split a 30-bit integer & transform
var u2cw = function(x){
    var p = new Array((x>>25)&31 , (x>>20)&31 , (x>>15)&31 , (x>>10)&31 , (x>>5)&31 , x&31);
    p.push(wsum1(p));
    p.push(wsum2(p));
    return p
};

// PIDgen
var PIDgen = function(x){
    x = x & NN_1;
    var p = new Array(8);
    var y = rndext(x);
    var z = encr(y);
    var cw = u2cw(z);
    for (j=0; j<=7; j++) p[j] = sigma[cw[j]];
    return p.join("");
};

var PIDArray = function(q){
    var out = new Array();
    for (let j=1; j<=q; j++){
        out.push(PIDgen(j));
    };
    return out;
};