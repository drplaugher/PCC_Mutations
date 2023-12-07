numvars=69;
--variables x and controls up=u+, um=u-
X=apply(numvars,i->"x"|(i+1));
UP=apply(numvars,i->"up"|(i+1));
UM=apply(numvars,i->"um"|(i+1));

--define ring
DEN=apply(join(X,UP,UM),v->v|"^2-"|v);
R=ZZ/2[join(X,UP,UM)/value,MonomialOrder=>Lex]/ideal(DEN/value);
X=X/value; UP=UP/value; UM=UM/value; 


RingElement | RingElement :=(x,y)->x+y+x*y;
RingElement & RingElement :=(x,y)->x*y;

-- Cytokines 
f1 = x22 | x45 | x50;
f2 = x55;
f3 = x29 | x53;
f4 = x4;
f5 = x19 | x44;
f6 = x6
f7 = x19 | x44
f8 = x8

--Stellate Cell (Blue)
f9 = x1
f10 = x2
f11 = x3
f12 = x4
f13 = x5
f14 = x6
f15 = x7
f16 = x8
f17 = x9 | x10 | x11
f18 = x13 | x17
f19 = x15
f20 = x16
f21 = x17
f22 = x12 | x28
f23 = x24 | x26
f24 = x21
f25 = (x27+1) & x18
f26 = (x32+1) & x23
f27 = x26
f28 = x25
f29 = x13 | x24
f30 = x29
f31 = x26
f32 = x28 | x26
f33 = x26
f34 = (x31+1) & x30
f35 = x28 & x30 & x36
f36 = x19 & ((x20 +1) | (x14+1)) & (x22 | x30)

--Pancreatic Cell (cyan)
f37 = x2 | x40
f38 = x3
f39 = x7
f40 = x40
f41 = x40
f42 = x37 | x43

-- RAS
f43 = x37 | x38
--f43 = x43+x43+1 -- mutated

-- SMAD
f44 = x39
--f44 = x44+x44 -- mutated

f45 = x41
f46 = (x52+1) & x42
f47 = x43
f48 = (x45+1) & (x44 | x64)
f49 = x47
f50 = x51
f51 = x46
f52 = x64
f53 = x49
f54 = (x57+1)
f55 = x53 | x59

-- CyclinD
f56 = (x48+1) & x50
--f56=x56+x56 -- mutated

f57 = (x56+1)
f58 = (x64+1) & (x50 | x45 | x51 | x59)
f59 = x49
f60 = (x55+1) & x51
f61 = (x58+1)
f62 = (x58+1) & (x66+1)
f63 = (x54+1) & (x51 | x64)

-- P53p 
f64 = (x63+1)
--f64 = x64+x64 --mutated

f65 = (x48+1) & x54
f66 = (x50+1) & (x64 | x62 | x61)
f67 = (x60+1) & (x50 | x62)
f68 = x66
f69 = x65 & (x59 | x55)

--define Boolean network
F = { f1 , f2 , f3 , f4 , f5 , f6 , f7 , f8 , f9 , f10 , f11 , f12 , f13 , f14 , f15 ,
 f16 , f17 , f18 , f19 , f20 , f21 , f22 , f23 , f24 , f25 , f26 , f27 , f28 , f29 ,
 f30 , f31 , f32 , f33 , f34 , f35 , f36 , f37 , f38 , f39 , f40 , f41 , f42 , f43 ,
 f44 , f45 , f46 , f47 , f48 , f49 , f50 , f51 , f52 , f53 , f54 , f55 , f56 , f57 ,
 f58 , f59 , f60 , f61 , f62 , f63 , f64 , f65 , f66 , f67 , f68 , f69 }; 



--define BN with control parameters
Fcontrol=apply(numvars,k->(UP_k+UM_k+1)*F_k+UP_k);

-- create equations
--excludig the cytokines TP53-64(off)	CDKN2A-56(off)	SMAD4-44(off)	KRAS-43(on), set restrictions
I=ideal( join(Fcontrol-X, {x43+1, x64, x56,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x57,x58,x59,x60,x61,x62,x63,x65,x66+1,x67,x68+1,x69} ));

GB=gens gb I
GB1=flatten entries GB

--Find generators that only contain parameters
EQUATIONS=select(GB1,P->(S:=support(P); set {}===(set S)*set(X) and P!=0));
Ring2=ZZ/2[join(X,UP,UM)/value];
CONTROLS=(EQUATIONS/toString/value)/factor;

--list all control policies
toString toList CONTROLS









######################################## ONLY UP ###############################################################
numvars=69;
--variables x and controls up=u+, um=u-
X=apply(numvars,i->"x"|(i+1));
UP=apply(numvars,i->"up"|(i+1));


--define ring
DEN=apply(join(X,UP),v->v|"^2-"|v);
R=ZZ/2[join(X,UP)/value,MonomialOrder=>Lex]/ideal(DEN/value);
X=X/value; UP=UP/value; 


RingElement | RingElement :=(x,y)->x+y+x*y;
RingElement & RingElement :=(x,y)->x*y;


--define Boolean network using functions above
F = { f1 , f2 , f3 , f4 , f5 , f6 , f7 , f8 , f9 , f10 , f11 , f12 , f13 , f14 , f15 ,
 f16 , f17 , f18 , f19 , f20 , f21 , f22 , f23 , f24 , f25 , f26 , f27 , f28 , f29 ,
 f30 , f31 , f32 , f33 , f34 , f35 , f36 , f37 , f38 , f39 , f40 , f41 , f42 , f43 ,
 f44 , f45 , f46 , f47 , f48 , f49 , f50 , f51 , f52 , f53 , f54 , f55 , f56 , f57 ,
 f58 , f59 , f60 , f61 , f62 , f63 , f64 , f65 , f66 , f67 , f68 , f69 }; 

--define BN with control parameters
Fcontrol=apply(numvars,k->(UP_k+1)*F_k+UP_k);


-- create equations
--excludig the cytokines TP53-64(off)	CDKN2A-56(off)	SMAD4-44(off)	KRAS-43(on), set restrictions
I=ideal( join(Fcontrol-X,{x43 +1, x64, x9 , x10 , x11 , x12 , x13 , x14 , x15 , x16 , x17 , x18 , x19 , x20 , x21 , x22 , x23 , x24 , x25 , x26+1 , x27 , x28 , x29 , x30 , x31 , x32 , x33+1 , x34 , x35 , x36 , x37 , x38 , x39 , x40 , x41 , x42 , x45 , x46 , x47, x49 , x50 , x51 , x53 , x54 , x55 , x57 , x58 , x59 , x60 , x61+1 , x62 , x65 , x66+1 , x67 , x68+1 , x69 }));



GB=gens gb I
GB1=flatten entries GB

--Find generators that only contain parameters
EQUATIONS=select(GB1,P->(S:=support(P); set {}===(set S)*set(X) and P!=0));
Ring2=ZZ/2[join(X,UP)/value];
CONTROLS=(EQUATIONS/toString/value)/factor;

--list all control policies
toString toList CONTROLS



	  
####################################### ONLY UM########################################################	  
numvars=69;
--variables x and controls up=u+, um=u-
X=apply(numvars,i->"x"|(i+1));
UM=apply(numvars,i->"um"|(i+1));

--define ring
DEN=apply(join(X,UM),v->v|"^2-"|v);
R=ZZ/2[join(X,UM)/value,MonomialOrder=>Lex]/ideal(DEN/value);
X=X/value; UM=UM/value; 


RingElement | RingElement :=(x,y)->x+y+x*y;
RingElement & RingElement :=(x,y)->x*y;


--define Boolean network using functions above
F = { f1 , f2 , f3 , f4 , f5 , f6 , f7 , f8 , f9 , f10 , f11 , f12 , f13 , f14 , f15 ,
 f16 , f17 , f18 , f19 , f20 , f21 , f22 , f23 , f24 , f25 , f26 , f27 , f28 , f29 ,
 f30 , f31 , f32 , f33 , f34 , f35 , f36 , f37 , f38 , f39 , f40 , f41 , f42 , f43 ,
 f44 , f45 , f46 , f47 , f48 , f49 , f50 , f51 , f52 , f53 , f54 , f55 , f56 , f57 ,
 f58 , f59 , f60 , f61 , f62 , f63 , f64 , f65 , f66 , f67 , f68 , f69 }; 
 

--define BN with control parameters
Fcontrol=apply(numvars,k->(UM_k+1)*F_k);


-- create equations
--excludig the cytokines TP53-64(off)	CDKN2A-56(off)	SMAD4-44(off)	KRAS-43(on), set restrictions
I=ideal( join(Fcontrol-X,{x43 +1, x64, x9+1 , x10 +1, x11+1  , x17 , x18 , x21+1 , x22+1 , x23+1 , x24+1 , x25 +1, x26+1 , x27 , x28 , x29 , x30 , x31 , x32 , x33+1 , x34 , x35 , x36 , x37 , x38 , x39 , x40 , x41 , x42 , x45 , x46 , x47, x49 , x50 , x51 , x53 , x54 , x55 , x58 , x59 , x60 , x61+1 , x62 , x65 , x66+1 , x67 , x68+1 , x69 }));


GB=gens gb I
GB1=flatten entries GB

--Find generators that only contain parameters
EQUATIONS=select(GB1,P->(S:=support(P); set {}===(set S)*set(X) and P!=0));
Ring2=ZZ/2[join(X,UM)/value];
CONTROLS=(EQUATIONS/toString/value)/factor;

--list all control policies
toString toList CONTROLS