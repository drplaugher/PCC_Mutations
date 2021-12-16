n=69;
--for edge control
VARSstring={{"x22","x45","x50"},
{"x55"},
{"x29","x53"},
{"x4"},
{"x19","x44"},
{"x6"},
{"x19","x44"},
{"x8"},
{"x1"},
{"x2"},
{"x3"},
{"x4"},
{"x5"},
{"x6"},
{"x7"},
{"x8"},
{"x9","x10","x11"},
{"x13","x17"},
{"x15"},
{"x16"},
{"x17"},
{"x12","x28"},
{"x24","x26"},
{"x21"},
{"x27","x18"},
{"x32","x23"},
{"x26"},
{"x25"},
{"x13","x24"},
{"x29"},
{"x26"},
{"x28","x26"},
{"x26"},
{"x31","x30"},
{"x28","x30","x36"},
{"x19","x20","x14","x22","x30"},
{"x2","x40"},
{"x3"},
{"x7"},
{"x40"},
{"x40"},
{"x37","x43"},
{"x37","x38"},
{"x39"},
{"x41"},
{"x52","x42"},
{"x43"},
{"x45","x44","x64"},
{"x47"},
{"x51"},
{"x46"},
{"x64"},
{"x49"},
{"x57"},
{"x53","x59"},
{"x48","x50"},
{"x56"},
{"x64","x50","x45","x51","x59"},
{"x49"},
{"x55","x51"},
{"x58"},
{"x58","x66"},
{"x54","x51","x64"},
{"x63"},
{"x48","x54"},
{"x50","x64","x62","x61"},
{"x60","x50","x62"},
{"x66"},
{"x65","x59","x55"}
};


--uji denote parameters from xi to xj
-- (for edge deletions)
PARSstring={{"u_(22,1)","u_(45,1)","u_(50,1)"},
{"u_(55,2)"},
{"u_(29,3)","u_(53,3)"},
{"u_(4,4)"},
{"u_(19,5)","u_(44,5)"},
{"u_(6,6)"},
{"u_(19,7)","u_(44,7)"},
{"u_(8,8)"},
{"u_(1,9)"},
{"u_(2,10)"},
{"u_(3,11)"},
{"u_(4,12)"},
{"u_(5,13)"},
{"u_(6,14)"},
{"u_(7,15)"},
{"u_(8,16)"},
{"u_(9,17)","u_(10,17)","u_(11,17)"},
{"u_(13,18)","u_(17,18)"},
{"u_(15,19)"},
{"u_(16,20)"},
{"u_(17,21)"},
{"u_(12,22)","u_(28,22)"},
{"u_(24,23)","u_(26,23)"},
{"u_(21,24)"},
{"u_(27,25)","u_(18,25)"},
{"u_(32,26)","u_(23,26)"},
{"u_(26,27)"},
{"u_(25,28)"},
{"u_(13,29)","u_(24,29)"},
{"u_(29,30)"},
{"u_(26,31)"},
{"u_(28,32)","u_(26,32)"},
{"u_(26,33)"},
{"u_(31,34)","u_(30,34)"},
{"u_(28,35)","u_(30,35)","u_(36,35)"},
{"u_(19,36)","u_(20,36)","u_(14,36)","u_(22,36)","u_(30,36)"},
{"u_(2,37)","u_(40,37)"},
{"u_(3,38)"},
{"u_(7,39)"},
{"u_(40,40)"},
{"u_(40,41)"},
{"u_(37,42)","u_(43,42)"},
{"u_(37,43)","u_(38,43)"},
{"u_(39,44)"},
{"u_(41,45)"},
{"u_(52,46)","u_(42,46)"},
{"u_(43,47)"},
{"u_(45,48)","u_(44,48)","u_(64,48)"},
{"u_(47,49)"},
{"u_(51,50)"},
{"u_(46,51)"},
{"u_(64,52)"},
{"u_(49,53)"},
{"u_(57,54)"},
{"u_(53,55)","u_(59,55)"},
{"u_(48,56)","u_(50,56)"},
{"u_(56,57)"},
{"u_(64,58)","u_(50,58)","u_(45,58)","u_(51,58)","u_(59,58)"},
{"u_(49,59)"},
{"u_(55,60)","u_(51,60)"},
{"u_(58,61)"},
{"u_(58,62)","u_(66,62)"},
{"u_(54,63)","u_(51,63)","u_(64,63)"},
{"u_(63,64)"},
{"u_(48,65)","u_(54,65)"},
{"u_(50,66)","u_(64,66)","u_(62,66)","u_(61,66)"},
{"u_(60,67)","u_(50,67)","u_(62,67)"},
{"u_(66,68)"},
{"u_(65,69)","u_(59,69)","u_(55,69)"}
};

PARSselected=PARSstring;

-- Code to find controllers
Xstring=apply(n,i->"x"|(i+1));
DEN=apply(join(Xstring,flatten PARSstring),v->v|"^2-"|v)
R=ZZ/2[join(Xstring,flatten PARSstring)/value]/ideal(DEN/value);
PARS=apply(n, i-> (PARSselected_i)/value);
VARS=apply(n, i-> (VARSstring_i)/value);

RingElement | RingElement :=(x,y)->x+y+x*y;
RingElement & RingElement :=(x,y)->x*y;


-- TP53-64(off)	CDKN2A-56(off)	SMAD4-44(off)	KRAS-43(on)
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




--encoding the edge controls
H=F;
HPAR={};
for k from 0 to n-1 do ( X:=VARS_k; P:=PARS_k;
    SUBS:=apply(#X, i->(
      X_i=>(P_i+1)*X_i
--	X_i=>(P_i+1)*X_i+P_i
        )
    );
    HPAR=append(HPAR,sub(H_k,SUBS));
)

-- define good and bad states TP53-64(off)	CDKN2A-56(off)	SMAD4-44(off)	KRAS-43(on)
GOODSTATE={{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0} };
BADSTATE={{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,x44,1,1,1,1,1,1,1,1,1,1,1,x56,1,1,1,1,1,1,1,0,1,0,1,0,1} 
};

--for blocking
EQNS2={}; 
for k from 0 to #BADSTATE-1 do (x:=BADSTATE_k;
    S:=apply(n,i->(Xstring_i=>x_i));
    HPARx:=apply(HPAR,hi->sub(hi,S)+1);
    EQNS2=join(EQNS2,flatten(HPARx-x));    
); 
EQNS2
I2=ideal(EQNS2);
gens gb I2

--for creating a new fixed point
EQNS3={}; 
for k from 0 to #GOODSTATE-1 do (x:=GOODSTATE_k;
    S:=apply(n,i->(Xstring_i=>x_i));
    HPARx:=apply(HPAR,hi->sub(hi,S));
    EQNS3=join(EQNS3,flatten(HPARx-x));    
); 
EQNS3
I3=ideal(EQNS3);
gens gb I3


GB1=flatten entries GB
--Find generators that only contain parameters
EQUATIONS=select(GB1,P->(S:=support(P); set {}===(set S)*set(X) and P!=0));
Ring2=ZZ/2[join(X,UP,UM)/value];
CONTROLS=(EQUATIONS/toString/value)/factor;

--list all control policies
toString toList CONTROLS
