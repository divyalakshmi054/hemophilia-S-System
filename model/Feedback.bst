// s-law model for extrinsic and intrinsic coagulation cascade -

#static::start
    TFPI,AT,TF,TM,PC,TRAUMA,SURFACE,HEM
#static::end

#dynamic::start
    FII,FVII,FV,FX,FVIII,FIX,FXI,FXII,FIIa,FVIIa,FVa,FXa,FVIIIa,FIXa,FXIa,FXIIa,FIIa_inactive,AP,PL
#dynamic::end

#structure::start
    // name::{reactants} --> {products}
    // 1. name: reaction name (used in rate section)
    // 2. {reactants}: list of dynamic species acting as a reactants
    // 3. {products}: list of dynamic species acting as a reactants

    // reactions -

    r1::FVII --> FVIIa
    r2::FVIII --> FVIIIa
    r3::FIX --> FIXa
    r4::FV --> FVa
    r5::FX --> FXa
    r6::FII --> FIIa
    r7::FXI --> FXIa
    r8::FXII --> FXIIa
    r9::FIIa --> FIIa_inactive
    r10::PL --> AP
#structure::end

#rate::start
    // name::{factor_1,factor_2,....}
    // 1: name = reaction name (must be in the #pragma::structure section)
    // 2: list of factor symbols (can be both static and dynamic, order does not matter)
    // 3: the order of the records does not matter here, we establish the order from the structure section
    r1::{FVII,TRAUMA,TF,TFPI}
    r2::{FVIII,AP,FIIa,HEM}
    r3::{FIX,AP,FVIIa,FXIa,HEM}
    r4::{FV,AP,FIIa}
    r5::{FX,AP,FVIIa,HEM}
    r6::{FII,AP,FXa,FVa,HEM}
    r7::{FXI,FXIIa}
    r8::{FXII,SURFACE}
    r9::{FIIa,AT}
    r10::{PL,FIIa}
#rate::end
