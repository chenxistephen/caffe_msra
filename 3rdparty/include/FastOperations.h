#ifndef FAST_OPERATIONS_H
#define FAST_OPERATIONS_H

#pragma warning(push, 4)
#pragma warning(disable: 4731)

namespace ns_base
{
   namespace
   {
   #if !defined _WIN64
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a city-block distance between two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE double CityBlockDist(const double* v1, const double* v2, size_t size_in)
      {
         double ret_val;
         __asm
         {
            mov eax, size_in
            shl eax,3

            sub  esp,12
            mov  [esp+4] ,esi
            mov  [esp] ,edi
                  
            mov esi, v1
            mov edi, v2
            
            //Prepare mask for the fabs
            pcmpeqw xmm3, xmm3 // All 1's
            psrlq   xmm3, 1    // Shift out the highest bit
            
            pxor xmm2,xmm2
            cmp  eax,32
            jl   simple_variant

            mov  ecx, eax
            and  ecx, 0xFFFFFFE0
            sub  eax, ecx
            lea  esi,[esi+ecx]
            lea  edi,[edi+ecx]
            neg  ecx
            
            sub_loop:
               movapd      xmm0,xmmword ptr [esi+ecx]
               subpd       xmm0,xmmword ptr [edi+ecx]
               andpd       xmm0,xmm3 // fabs
               addpd       xmm2,xmm0

               movapd      xmm1,xmmword ptr [esi+ecx+16]
               subpd       xmm1,xmmword ptr [edi+ecx+16]
               andpd       xmm1,xmm3 // fabs
               addpd       xmm2,xmm1
               add         ecx,32
            jnz sub_loop
            
            test   eax,0x10
            jz     no_2
            movapd      xmm0,xmmword ptr [esi+ecx]
            subpd       xmm0,xmmword ptr [edi+ecx]
            andpd       xmm0,xmm3 // fabs
            addpd       xmm2,xmm0
            add         ecx,16
         
         no_2:
            test   eax,0x8
            jz     epilog

            movsd       xmm0,xmmword ptr [esi+ecx]
            subsd       xmm0,xmmword ptr [edi+ecx]
            andpd       xmm0,xmm3 // fabs
            addsd       xmm2,xmm0
            jmp         epilog

         simple_variant:
            xor    ecx,ecx
            test   eax,0x10
            jz     no_2_simple
            movapd      xmm0,xmmword ptr [esi+ecx]
            subpd       xmm0,xmmword ptr [edi+ecx]
            andpd       xmm0,xmm3 // fabs
            addpd       xmm2,xmm0
            add         ecx,16
         
         no_2_simple:
            test   eax,0x8
            jz     epilog

            movsd       xmm0,xmmword ptr [esi+ecx]
            subsd       xmm0,xmmword ptr [edi+ecx]
            andpd       xmm0,xmm3 // fabs
            addsd       xmm2,xmm0

         epilog:
            mov  esi,[esp+4]
            mov  edi,[esp]
            add  esp,12
            
            movhlps     xmm0,xmm2
            addsd       xmm2,xmm0
            movsd       ret_val,xmm2;
         }
         return ret_val;
      }
      
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates an euclidean distance between two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE double EuclideanDist(const double* v1, const double* v2, size_t size_in)
      {
         double ret_val;
         __asm
         {
            mov eax, size_in
            shl eax,3

            sub  esp,12
            mov  [esp+4] ,esi
            mov  [esp] ,edi
                  
            mov esi, v1
            mov edi, v2
            
            pxor xmm2,xmm2
            cmp  eax,32
            jl   simple_variant

            mov  ecx, eax
            and  ecx, 0xFFFFFFE0
            sub  eax, ecx
            lea  esi,[esi+ecx]
            lea  edi,[edi+ecx]
            neg  ecx
            
            sub_loop:
               movapd      xmm0,xmmword ptr [esi+ecx]
               subpd       xmm0,xmmword ptr [edi+ecx]
               mulpd       xmm0,xmm0
               addpd       xmm2,xmm0

               movapd      xmm1,xmmword ptr [esi+ecx+16]
               subpd       xmm1,xmmword ptr [edi+ecx+16]
               mulpd       xmm1,xmm1
               addpd       xmm2,xmm1
               add         ecx,32
            jnz sub_loop
            
            test   eax,0x10
            jz     no_2
            movapd      xmm0,xmmword ptr [esi+ecx]
            subpd       xmm0,xmmword ptr [edi+ecx]
            mulpd       xmm0,xmm0
            addpd       xmm2,xmm0
            add         ecx,16
         
         no_2:
            test   eax,0x8
            jz     epilog

            movsd       xmm0,xmmword ptr [esi+ecx]
            subsd       xmm0,xmmword ptr [edi+ecx]
            mulsd       xmm0,xmm0
            addsd       xmm2,xmm0
            jmp         epilog

         simple_variant:
            xor    ecx,ecx
            test   eax,0x10
            jz     no_2_simple
            movapd      xmm0,xmmword ptr [esi+ecx]
            subpd       xmm0,xmmword ptr [edi+ecx]
            mulpd       xmm0,xmm0
            addpd       xmm2,xmm0
            add         ecx,16
         
         no_2_simple:
            test   eax,0x8
            jz     epilog

            movsd       xmm0,xmmword ptr [esi+ecx]
            subsd       xmm0,xmmword ptr [edi+ecx]
            mulsd       xmm0,xmm0
            addsd       xmm2,xmm0

         epilog:
            mov  esi,[esp+4]
            mov  edi,[esp]
            add  esp,12
            
            movhlps     xmm0,xmm2
            addsd       xmm2,xmm0
            movsd       ret_val,xmm2;
         }
         return ret_val;
      }

      FORCE_INLINE float EuclideanDistFloat(const float* v1, const float* v2, size_t size_in)
      {
         double ret_val = 0;
         for (size_t i=0; i<size_in; ++i)
         {
            double d = v1[i]-v2[i];
            ret_val += d*d;
         }
         return (float)ret_val;
      }

      FORCE_INLINE void MulSquareMatrixByVector(const double* m, const double* v, size_t size_in, double* res)
      {
         __asm
         {
            mov eax, size_in
            shl eax,3

            sub  esp,16
            mov  [esp+12],ebx
            mov  [esp+8] ,esi
            mov  [esp+4] ,edi
            mov  [esp]   ,ebp
            
            mov esi, m
            mov edi, v
            mov ebp, res
            lea edi,[edi+eax]
            lea ebp,[ebp+eax]

            cmp  eax,64
            jl   simple_variant
            test eax,31
            jz   devisible_by_4
            test eax,15
            jz   devisible_by_2
            test eax,16
            jz   remainder_1
               
               neg eax
               lea ebx,[eax+8]
               rem3_outer_loop:
                  //first line
                  sub         esi,eax
                  lea         ecx,[eax+32+24]
                  movapd      xmm0,xmmword ptr [esi+eax]
                  mulpd       xmm0,xmmword ptr [edi+eax]
                  movapd      xmm1,xmmword ptr [esi+eax+16]
                  mulpd       xmm1,xmmword ptr [edi+eax+16]

                  rem3_inner_loop1:
                     movapd      xmm2,xmmword ptr [esi+ecx-24]
                     mulpd       xmm2,xmmword ptr [edi+ecx-24]
                     addpd       xmm0,xmm2
                     movapd      xmm3,xmmword ptr [esi+ecx-8]
                     mulpd       xmm3,xmmword ptr [edi+ecx-8]
                     addpd       xmm1,xmm3
                  add         ecx,32
                  jnz rem3_inner_loop1

                  movapd      xmm2,xmmword ptr [esi+ecx-24]
                  mulpd       xmm2,xmmword ptr [edi+ecx-24]
                  addpd       xmm0,xmm2

                  movsd       xmm3,xmmword ptr [esi+ecx-8]
                  mulsd       xmm3,xmmword ptr [edi+ecx-8]
                  addsd       xmm1,xmm3

                  addpd       xmm0,xmm1
                  movhlps     xmm1,xmm0
                  addsd       xmm0,xmm1
                  movsd       [ebp+ebx-8],xmm0

                  // second line
                  sub         esi,eax
                  lea         ecx,[eax+32+24]
                  movq        xmm0,qword ptr [esi+eax]
                  movhps      xmm0,qword ptr [esi+eax+8]
                  mulpd       xmm0,xmmword ptr [edi+eax]
                  movq        xmm1,qword ptr [esi+eax+16]
                  movhps      xmm1,qword ptr [esi+eax+24]
                  mulpd       xmm1,xmmword ptr [edi+eax+16]

                  rem3_inner_loop2:
                     movq        xmm2,qword ptr [esi+ecx-24]
                     movhps      xmm2,qword ptr [esi+ecx-16]
                     mulpd       xmm2,xmmword ptr [edi+ecx-24]
                     addpd       xmm0,xmm2
                     movq        xmm3,qword ptr [esi+ecx-8]
                     movhps      xmm3,qword ptr [esi+ecx]
                     mulpd       xmm3,xmmword ptr [edi+ecx-8]
                     addpd       xmm1,xmm3
                  add         ecx,32
                  jnz rem3_inner_loop2

                  movq        xmm2,qword ptr [esi+ecx-24]
                  movhps      xmm2,qword ptr [esi+ecx-16]
                  mulpd       xmm2,xmmword ptr [edi+ecx-24]
                  addpd       xmm0,xmm2

                  movsd       xmm3,xmmword ptr [esi+ecx-8]
                  mulsd       xmm3,xmmword ptr [edi+ecx-8]
                  addsd       xmm1,xmm3

                  addpd       xmm0,xmm1
                  movhlps     xmm1,xmm0
                  addsd       xmm0,xmm1
                  movsd       [ebp+ebx],xmm0
               add         ebx,16
               jnz rem3_outer_loop
               
               //last line
               sub         esi,eax
               lea         ecx,[eax+32+24]
               movapd      xmm0,xmmword ptr [esi+eax]
               mulpd       xmm0,xmmword ptr [edi+eax]
               movapd      xmm1,xmmword ptr [esi+eax+16]
               mulpd       xmm1,xmmword ptr [edi+eax+16]

               rem3_inner_loop_last:
                  movapd      xmm2,xmmword ptr [esi+ecx-24]
                  mulpd       xmm2,xmmword ptr [edi+ecx-24]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx-8]
                  mulpd       xmm3,xmmword ptr [edi+ecx-8]
                  addpd       xmm1,xmm3
               add         ecx,32
               jnz rem3_inner_loop_last

               movapd      xmm2,xmmword ptr [esi+ecx-24]
               mulpd       xmm2,xmmword ptr [edi+ecx-24]
               addpd       xmm0,xmm2

               movsd       xmm3,xmmword ptr [esi+ecx-8]
               mulsd       xmm3,xmmword ptr [edi+ecx-8]
               addsd       xmm1,xmm3

               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+ebx-8],xmm0
            jmp epilog

            devisible_by_2:
               neg eax
               mov ebx,eax
               d2_outer_loop:
                  sub         esi,eax
                  lea         ecx,[eax+32+16]
                  movapd      xmm0,xmmword ptr [esi+eax]
                  mulpd       xmm0,xmmword ptr [edi+eax]
                  movapd      xmm1,xmmword ptr [esi+eax+16]
                  mulpd       xmm1,xmmword ptr [edi+eax+16]
                  
                  d2_inner_loop:
                     movapd      xmm2,xmmword ptr [esi+ecx-16]
                     mulpd       xmm2,xmmword ptr [edi+ecx-16]
                     addpd       xmm0,xmm2
                     movapd      xmm3,xmmword ptr [esi+ecx]
                     mulpd       xmm3,xmmword ptr [edi+ecx]
                     addpd       xmm1,xmm3
                  add         ecx,32
                  jnz d2_inner_loop

                  movapd      xmm2,xmmword ptr [esi+ecx-16]
                  mulpd       xmm2,xmmword ptr [edi+ecx-16]
                  addpd       xmm0,xmm2

                  addpd       xmm0,xmm1
                  movhlps     xmm1,xmm0
                  addsd       xmm0,xmm1
                  movsd       [ebp+ebx],xmm0
               add         ebx,8
               jnz d2_outer_loop
            jmp epilog

            remainder_1:
               neg eax
               lea ebx,[eax+8]
               rem1_outer_loop:
                  //first line
                  sub         esi,eax
                  lea         ecx,[eax+32+8]
                  movapd      xmm0,xmmword ptr [esi+eax]
                  mulpd       xmm0,xmmword ptr [edi+eax]
                  movapd      xmm1,xmmword ptr [esi+eax+16]
                  mulpd       xmm1,xmmword ptr [edi+eax+16]

                  rem1_inner_loop1:
                     movapd      xmm2,xmmword ptr [esi+ecx-8]
                     mulpd       xmm2,xmmword ptr [edi+ecx-8]
                     addpd       xmm0,xmm2
                     movapd      xmm3,xmmword ptr [esi+ecx+8]
                     mulpd       xmm3,xmmword ptr [edi+ecx+8]
                     addpd       xmm1,xmm3
                  add         ecx,32
                  jnz rem1_inner_loop1

                  movsd       xmm2,xmmword ptr [esi+ecx-8]
                  mulsd       xmm2,xmmword ptr [edi+ecx-8]
                  addsd       xmm0,xmm2

                  addpd       xmm0,xmm1
                  movhlps     xmm1,xmm0
                  addsd       xmm0,xmm1
                  movsd       [ebp+ebx-8],xmm0

                  // second line
                  sub         esi,eax
                  lea         ecx,[eax+32+8]
                  movq        xmm0,qword ptr [esi+eax]
                  movhps      xmm0,qword ptr [esi+eax+8]
                  mulpd       xmm0,xmmword ptr [edi+eax]
                  movq        xmm1,qword ptr [esi+eax+16]
                  movhps      xmm1,qword ptr [esi+eax+24]
                  mulpd       xmm1,xmmword ptr [edi+eax+16]

                  rem1_inner_loop2:
                     movq        xmm2,qword ptr [esi+ecx-8]
                     movhps      xmm2,qword ptr [esi+ecx]
                     mulpd       xmm2,xmmword ptr [edi+ecx-8]
                     addpd       xmm0,xmm2
                     movq        xmm3,qword ptr [esi+ecx+8]
                     movhps      xmm3,qword ptr [esi+ecx+16]
                     mulpd       xmm3,xmmword ptr [edi+ecx+8]
                     addpd       xmm1,xmm3
                  add         ecx,32
                  jnz rem1_inner_loop2

                  movsd       xmm2,xmmword ptr [esi+ecx-8]
                  mulsd       xmm2,xmmword ptr [edi+ecx-8]
                  addsd       xmm0,xmm2

                  addpd       xmm0,xmm1
                  movhlps     xmm1,xmm0
                  addsd       xmm0,xmm1
                  movsd       [ebp+ebx],xmm0

               add         ebx,16
               jnz rem1_outer_loop
               
               //last line
               sub         esi,eax
               lea         ecx,[eax+32+8]
               movapd      xmm0,xmmword ptr [esi+eax]
               mulpd       xmm0,xmmword ptr [edi+eax]
               movapd      xmm1,xmmword ptr [esi+eax+16]
               mulpd       xmm1,xmmword ptr [edi+eax+16]

               rem1_inner_loop_last:
                  movapd      xmm2,xmmword ptr [esi+ecx-8]
                  mulpd       xmm2,xmmword ptr [edi+ecx-8]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+8]
                  mulpd       xmm3,xmmword ptr [edi+ecx+8]
                  addpd       xmm1,xmm3
               add         ecx,32
               jnz rem1_inner_loop_last

               movsd       xmm2,xmmword ptr [esi+ecx-8]
               mulsd       xmm2,xmmword ptr [edi+ecx-8]
               addsd       xmm0,xmm2

               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+ebx-8],xmm0
            jmp epilog

            devisible_by_4:
               neg eax
               mov ebx,eax

               d4_outer_loop:
                  sub         esi,eax
                  lea         ecx,[eax+32]
                  movapd      xmm0,xmmword ptr [esi+eax]
                  mulpd       xmm0,xmmword ptr [edi+eax]
                  movapd      xmm1,xmmword ptr [esi+eax+16]
                  mulpd       xmm1,xmmword ptr [edi+eax+16]
                  d4_inner_loop:
                     movapd      xmm2,xmmword ptr [esi+ecx]
                     mulpd       xmm2,xmmword ptr [edi+ecx]
                     addpd       xmm0,xmm2
                     movapd      xmm3,xmmword ptr [esi+ecx+16]
                     mulpd       xmm3,xmmword ptr [edi+ecx+16]
                     addpd       xmm1,xmm3
                  add         ecx,32
                  jnz d4_inner_loop

                  addpd       xmm0,xmm1
                  movhlps     xmm1,xmm0
                  addsd       xmm0,xmm1
                  movsd       [ebp+ebx],xmm0
               add         ebx,8
               jnz d4_outer_loop
            jmp epilog

         simple_variant:
               neg eax
               jz  epilog
               mov ebx,eax
               outer_loop_sv:
                  sub esi,eax
                  mov ecx,eax
                  fldz
                  inner_loop_sv:
                     fld   qword ptr[esi+ecx]
                     fmul  qword ptr[edi+ecx]
                     faddp st(1),st(0)
                     add ecx,8
                  jnz inner_loop_sv
                  fstp  qword ptr[ebp+ebx]
                  add ebx,8
               jnz outer_loop_sv
         epilog:
            
            mov  ebx,[esp+12]
            mov  esi,[esp+8]
            mov  edi,[esp+4]
            mov  ebp,[esp]
            add  esp,16
         }
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies rectangular matrix by a vector
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MulMatrixByVector(const double* m, const double* v, size_t n_rows, size_t n_cols, double* res)
      {
         __asm
         {
            mov eax, n_cols
            shl eax,3

            sub  esp,20
            mov  [esp+16],edx
            mov  [esp+12],ebx
            mov  [esp+8] ,esi
            mov  [esp+4] ,edi
            mov  [esp]   ,ebp
            
            mov esi, m
            mov edi, v
            mov ebx, n_rows
            mov ebp, res

            cmp  eax,64
            jl   simple_variant

            mov edx,eax
            and edx,0xFFFFFFE0
            lea edi,[edi+edx]
            lea ebp,[ebp+8*ebx]

            sub eax,edx
            neg edx
            neg ebx

            test eax,0x10
            jz   no_2

            test eax,0x8
            jz   outer_loop_even

            outer_loop:
               sub         esi,edx
               lea         ecx,[edx+32]
               movapd      xmm0,xmmword ptr [esi+edx]
               mulpd       xmm0,xmmword ptr [edi+edx]
               movapd      xmm1,xmmword ptr [esi+edx+16]
               mulpd       xmm1,xmmword ptr [edi+edx+16]
               inner_loop:
                  movapd      xmm2,xmmword ptr [esi+ecx]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+16]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz inner_loop
               
               movapd      xmm2,xmmword ptr [esi+ecx]
               mulpd       xmm2,xmmword ptr [edi+ecx]
               addpd       xmm0,xmm2
               add         ecx,16
               movsd       xmm3,xmmword ptr [esi+ecx]
               mulsd       xmm3,xmmword ptr [edi+ecx]
               addpd       xmm1,xmm3
               
               add         esi,eax
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+8*ebx],xmm0
               add         ebx,1
               jz          epilog

               //------------ unaligned part ----------------
               sub         esi,edx
               lea         ecx,[edx+32]

               movq        xmm0,qword ptr [esi+edx]
               movhps      xmm0,qword ptr [esi+edx+8]
               mulpd       xmm0,xmmword ptr [edi+edx]

               movq        xmm1,qword ptr [esi+edx+16]
               movhps      xmm1,qword ptr [esi+edx+24]
               mulpd       xmm1,xmmword ptr [edi+edx+16]

               inner_loop_u:
                  movq        xmm2,qword ptr [esi+ecx]
                  movhps      xmm2,qword ptr [esi+ecx+8]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2

                  movq        xmm3,qword ptr [esi+ecx+16]
                  movhps      xmm3,qword ptr [esi+ecx+24]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz inner_loop_u
               
               movq        xmm2,qword ptr [esi+ecx]
               movhps      xmm2,qword ptr [esi+ecx+8]
               mulpd       xmm2,xmmword ptr [edi+ecx]
               addpd       xmm0,xmm2
               add         ecx,16
               movsd       xmm3,xmmword ptr [esi+ecx]
               mulsd       xmm3,xmmword ptr [edi+ecx]
               addpd       xmm1,xmm3

               add         esi,eax
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+8*ebx],xmm0
               add         ebx,1
            jnz outer_loop
            jmp epilog

            // ============  even variant ================ 
            outer_loop_even:
               sub         esi,edx
               lea         ecx,[edx+32]
               movapd      xmm0,xmmword ptr [esi+edx]
               mulpd       xmm0,xmmword ptr [edi+edx]
               movapd      xmm1,xmmword ptr [esi+edx+16]
               mulpd       xmm1,xmmword ptr [edi+edx+16]
               inner_loop_even:
                  movapd      xmm2,xmmword ptr [esi+ecx]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+16]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz inner_loop_even
               
               movapd      xmm2,xmmword ptr [esi+ecx]
               mulpd       xmm2,xmmword ptr [edi+ecx]
               addpd       xmm0,xmm2
               add         esi,eax
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+8*ebx],xmm0
               add         ebx,1
            jnz outer_loop_even
            jmp epilog

            //========== No 2 Part ============

         no_2:
            test eax,0x8
            jz   outer_loop_n2_even

            outer_loop_n2:
               sub         esi,edx
               lea         ecx,[edx+32]
               movapd      xmm0,xmmword ptr [esi+edx]
               mulpd       xmm0,xmmword ptr [edi+edx]
               movapd      xmm1,xmmword ptr [esi+edx+16]
               mulpd       xmm1,xmmword ptr [edi+edx+16]
               inner_loop_n2:
                  movapd      xmm2,xmmword ptr [esi+ecx]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+16]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz inner_loop_n2
               
               movsd       xmm3,xmmword ptr [esi+ecx]
               mulsd       xmm3,xmmword ptr [edi+ecx]
               addpd       xmm1,xmm3
               
               add         esi,eax
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+8*ebx],xmm0
               add         ebx,1
               jz          epilog

               //------------ unaligned part ----------------
               sub         esi,edx
               lea         ecx,[edx+32]

               movq        xmm0,qword ptr [esi+edx]
               movhps      xmm0,qword ptr [esi+edx+8]
               mulpd       xmm0,xmmword ptr [edi+edx]

               movq        xmm1,qword ptr [esi+edx+16]
               movhps      xmm1,qword ptr [esi+edx+24]
               mulpd       xmm1,xmmword ptr [edi+edx+16]

               inner_loop_u_n2:
                  movq        xmm2,qword ptr [esi+ecx]
                  movhps      xmm2,qword ptr [esi+ecx+8]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2

                  movq        xmm3,qword ptr [esi+ecx+16]
                  movhps      xmm3,qword ptr [esi+ecx+24]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz inner_loop_u_n2
               movsd       xmm3,xmmword ptr [esi+ecx]
               mulsd       xmm3,xmmword ptr [edi+ecx]
               addpd       xmm1,xmm3

               add         esi,eax
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+8*ebx],xmm0
               add         ebx,1
            jnz outer_loop_n2
            jmp epilog

            //============== Even N2 Part ================

            outer_loop_n2_even:
               sub         esi,edx
               lea         ecx,[edx+32]
               movapd      xmm0,xmmword ptr [esi+edx]
               mulpd       xmm0,xmmword ptr [edi+edx]
               movapd      xmm1,xmmword ptr [esi+edx+16]
               mulpd       xmm1,xmmword ptr [edi+edx+16]
               inner_loop_n2_even:
                  movapd      xmm2,xmmword ptr [esi+ecx]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+16]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz inner_loop_n2_even
               add         esi,eax
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+8*ebx],xmm0
               add         ebx,1
            jnz outer_loop_n2_even
            jmp epilog

            //============================================

         simple_variant:
               lea edi,[edi+eax]
               lea ebp,[ebp+8*ebx]
               neg eax
               neg ebx
               jz  epilog         
               outer_loop_sv:
                  sub esi,eax
                  mov ecx,eax
                  fldz
                  inner_loop_sv:
                     fld   qword ptr[esi+ecx]
                     fmul  qword ptr[edi+ecx]
                     faddp st(1),st(0)
                     add ecx,8
                  jnz inner_loop_sv
                  fstp  qword ptr[ebp+8*ebx]
                  add ebx,1
               jnz outer_loop_sv
         epilog:
            mov  edx,[esp+16]
            mov  ebx,[esp+12]
            mov  esi,[esp+8]
            mov  edi,[esp+4]
            mov  ebp,[esp]
            add  esp,20
         }
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies rectangular matrix by a vector
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MulMatrixByVectorDF(const double* m, const double* v, size_t n_rows, size_t n_cols, float* res)
      {
         double sum;
         for (size_t i=0; i<n_rows; ++i)
         {
            sum = 0;
            for (size_t j=0; j<n_cols; ++j)
               sum += m[i*n_cols+j]*v[j];
            res[i] = (float)sum;
         }
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies rectangular matrix by a vector
      // Vector dimension should be divisible by 8
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MulMatrixByVectorFloat(const float* m, const float* v, size_t n_rows, size_t n_cols, float* res)
      {
         double sum;
         for (size_t i=0; i<n_rows; ++i)
         {
            sum = 0;
            for (size_t j=0; j<n_cols; ++j)
               sum += m[i*n_cols+j]*v[j];
            res[i] = (float)sum;
         }
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function is multiplying upper-triangular matrix by vector. Matrix and vector should be padded with 0 to the multiple of 4
      // Size_in is PADDED SIZE - not the original one
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MulTriangularMatrixByVector4(const double* m, const double* v, size_t size_in, double* res)
      {
         __asm
         {
            mov eax, size_in
            shl eax,3

            sub  esp,20
            mov  [esp+16],edx
            mov  [esp+12],ebx
            mov  [esp+8] ,esi
            mov  [esp+4] ,edi
            mov  [esp]   ,ebp
            
            mov esi, m
            mov edi, v
            mov ebp, res
            
            cmp  eax,64
            jl   simple_variant

            lea edi,[edi+eax]
            lea ebp,[ebp+eax-32]
            xor edx,edx
            neg eax
            lea ebx,[eax+32]

            d4_outer_loop:
               sub         esi,eax
               lea         ecx,[eax+edx+32]
               movapd      xmm0,xmmword ptr [esi+ecx-32]
               mulpd       xmm0,xmmword ptr [edi+ecx-32]
               movapd      xmm1,xmmword ptr [esi+ecx-16]
               mulpd       xmm1,xmmword ptr [edi+ecx-16]
               d4_inner_loop1:
                  movapd      xmm2,xmmword ptr [esi+ecx]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+16]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz d4_inner_loop1
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+ebx],xmm0
               add         ebx,8
               //------------------------------------------------
               sub         esi,eax
               lea         ecx,[eax+edx+32]
               movsd       xmm0,xmmword ptr [esi+ecx-24]
               mulsd       xmm0,xmmword ptr [edi+ecx-24]
               movapd      xmm1,xmmword ptr [esi+ecx-16]
               mulpd       xmm1,xmmword ptr [edi+ecx-16]
               d4_inner_loop2:
                  movapd      xmm2,xmmword ptr [esi+ecx]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+16]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz d4_inner_loop2
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+ebx],xmm0
               add         ebx,8
               //------------------------------------------------
               sub         esi,eax
               lea         ecx,[eax+edx+32]
               pxor        xmm0,xmm0
               movapd      xmm1,xmmword ptr [esi+ecx-16]
               mulpd       xmm1,xmmword ptr [edi+ecx-16]
               d4_inner_loop3:
                  movapd      xmm2,xmmword ptr [esi+ecx]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+16]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz d4_inner_loop3
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+ebx],xmm0
               add         ebx,8
               //------------------------------------------------
               sub         esi,eax
               lea         ecx,[eax+edx+32]
               pxor        xmm0,xmm0
               movsd       xmm1,xmmword ptr [esi+ecx-8]
               mulsd       xmm1,xmmword ptr [edi+ecx-8]
               d4_inner_loop4:
                  movapd      xmm2,xmmword ptr [esi+ecx]
                  mulpd       xmm2,xmmword ptr [edi+ecx]
                  addpd       xmm0,xmm2
                  movapd      xmm3,xmmword ptr [esi+ecx+16]
                  mulpd       xmm3,xmmword ptr [edi+ecx+16]
                  addpd       xmm1,xmm3
                  add         ecx,32
               jnz d4_inner_loop4
               addpd       xmm0,xmm1
               movhlps     xmm1,xmm0
               addsd       xmm0,xmm1
               movsd       [ebp+ebx],xmm0
               add         edx,32
               add         ebx,8
            jnz d4_outer_loop

            sub         esi,eax
            lea         ecx,[eax+edx]
            movapd      xmm0,xmmword ptr [esi+ecx]
            mulpd       xmm0,xmmword ptr [edi+ecx]
            movapd      xmm1,xmmword ptr [esi+ecx+16]
            mulpd       xmm1,xmmword ptr [edi+ecx+16]
            addpd       xmm0,xmm1
            movhlps     xmm1,xmm0
            addsd       xmm0,xmm1
            movsd       [ebp+ebx],xmm0
            //------------------------------------------------
            sub         esi,eax
            movsd       xmm0,xmmword ptr [esi+ecx+8]
            mulsd       xmm0,xmmword ptr [edi+ecx+8]
            movapd      xmm1,xmmword ptr [esi+ecx+16]
            mulpd       xmm1,xmmword ptr [edi+ecx+16]
            addpd       xmm0,xmm1
            movhlps     xmm1,xmm0
            addsd       xmm0,xmm1
            movsd       [ebp+ebx+8],xmm0
            //------------------------------------------------
            sub         esi,eax
            movapd      xmm0,xmmword ptr [esi+ecx+16]
            mulpd       xmm0,xmmword ptr [edi+ecx+16]
            movhlps     xmm1,xmm0
            addsd       xmm0,xmm1
            movsd       [ebp+ebx+16],xmm0
            //------------------------------------------------
            sub         esi,eax
            movsd       xmm0,xmmword ptr [esi+ecx+24]
            mulsd       xmm0,xmmword ptr [edi+ecx+24]
            movsd       [ebp+ebx+24],xmm0
            jmp epilog

         simple_variant:
            neg eax
            jz  epilog
            movapd      xmm0,xmmword ptr [esi]
            mulpd       xmm0,xmmword ptr [edi]
            movapd      xmm1,xmmword ptr [esi+16]
            mulpd       xmm1,xmmword ptr [edi+16]
            addpd       xmm0,xmm1
            movhlps     xmm1,xmm0
            addsd       xmm0,xmm1
            movsd       [ebp],xmm0

            movsd       xmm2,xmmword ptr [esi+40]
            mulsd       xmm2,xmmword ptr [edi+8]
            movapd      xmm3,xmmword ptr [esi+48]
            mulpd       xmm3,xmmword ptr [edi+16]
            addpd       xmm2,xmm3
            movhlps     xmm3,xmm2
            addsd       xmm2,xmm3
            movsd       [ebp+8],xmm2

            movapd      xmm4,xmmword ptr [esi+80]
            mulpd       xmm4,xmmword ptr [edi+16]
            movhlps     xmm5,xmm4
            addsd       xmm4,xmm5
            movsd       [ebp+16],xmm4

            movsd       xmm6,xmmword ptr [esi+120]
            mulsd       xmm6,xmmword ptr [edi+24]
            movsd       [ebp+24],xmm6
               
         epilog:
            mov  edx,[esp+16]
            mov  ebx,[esp+12]
            mov  esi,[esp+8]
            mov  edi,[esp+4]
            mov  ebp,[esp]
            add  esp,20
         }
      }
      
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function subtracts two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void SubVectors(const double* v1, const double* v2, size_t size_in, double* res)
      {
         __asm
         {
            mov eax, size_in
            shl eax,3

            sub  esp,12
            mov  [esp+8] ,esi
            mov  [esp+4] ,edi
            mov  [esp]   ,ebp
            
            mov esi, v1
            mov edi, v2
            mov ebp, res

            cmp  eax,32
            jl   simple_variant

            mov  ecx, eax
            and  ecx, 0xFFFFFFE0
            sub  eax, ecx
            lea  esi,[esi+ecx]
            lea  edi,[edi+ecx]
            lea  ebp,[ebp+ecx]
            neg  ecx

            sub_loop:
               movapd      xmm0,xmmword ptr [esi+ecx]
               subpd       xmm0,xmmword ptr [edi+ecx]
               movapd      xmmword ptr [ebp+ecx],xmm0

               movapd      xmm1,xmmword ptr [esi+ecx+16]
               subpd       xmm1,xmmword ptr [edi+ecx+16]
               movapd      xmmword ptr [ebp+ecx+16],xmm1
               add         ecx,32
            jnz sub_loop
            
            test   eax,0x10
            jz     no_2
            movapd      xmm0,xmmword ptr [esi+ecx]
            subpd       xmm0,xmmword ptr [edi+ecx]
            movapd      xmmword ptr [ebp+ecx],xmm0
            add         ecx,16
         
         no_2:
            test   eax,0x8
            jz     epilog

            movsd       xmm0,xmmword ptr [esi+ecx]
            subsd       xmm0,xmmword ptr [edi+ecx]
            movsd       xmmword ptr [ebp+ecx],xmm0
            jmp         epilog

         simple_variant:
            xor    ecx,ecx
            test   eax,0x10
            jz     no_2_simple
            movapd      xmm0,xmmword ptr [esi+ecx]
            subpd       xmm0,xmmword ptr [edi+ecx]
            movapd      xmmword ptr [ebp+ecx],xmm0
            add         ecx,16
         
         no_2_simple:
            test   eax,0x8
            jz     epilog

            movsd       xmm0,xmmword ptr [esi+ecx]
            subsd       xmm0,xmmword ptr [edi+ecx]
            movsd       xmmword ptr [ebp+ecx],xmm0

         epilog:
            mov  esi,[esp+8]
            mov  edi,[esp+4]
            mov  ebp,[esp]
            add  esp,12
         }
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function subtracts two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void SubVectorsU(const double* v1, const double* v2, size_t size_in, double* res)
      {
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] - v2[i];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates element wise product of two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultVectors(const double* v1, const double* v2, size_t size_in, double* res)
      {
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] * v2[i];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void AddVectorsIP(double* v1, const double* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i];
      }

      FORCE_INLINE void AddVectorsIP3(double* v1, const double* v2, const double* v3, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] + v3[i];
      }

      FORCE_INLINE void AddVectorsIP4(double* v1, const double* v2, const double* v3, const double* v4, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] + v3[i] + v4[i];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultByConstAndAdd(const double* v, double c, size_t size_in, double* res)
      {
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and places result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultByConst(const double* v, double c, size_t size_in, double* res)
      {
         for (size_t l=0; l<size_in; ++l)
            res[l] = c*v[l];
      }

      FORCE_INLINE void MultByConstIP(double* v, double c, size_t size_in)
      {
         for (size_t l=0; l<size_in; ++l)
            v[l] *= c;
      }

      FORCE_INLINE void AddConstIP(double* v, double c, size_t size_in)
      {
         for (size_t l=0; l<size_in; ++l)
            v[l] += c;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // size must be devisible by 4 and everything should be aligned to 16 bytes!!!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultByConstAndAddFloat4(const float* v, float c, size_t size_in, float* res)
      {
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void AddVectorsIPFloat(float* v1, const float* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i];
      }

      FORCE_INLINE void AddVectorsIP3Float(float* v1, const float* v2, const float* v3, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] + v3[i];
      }

      FORCE_INLINE void AddVectorsIP4Float(float* v1, const float* v2, const float* v3, const float* v4, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] + v3[i] + v4[i];
      }

      FORCE_INLINE void AddVectorSqr(double* v1, const double* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] * v2[i];
      }

      FORCE_INLINE void MaxVectorsIP(double* v1, const double* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::max(v1[i], v2[i]);
      }

      FORCE_INLINE void MinVectorsIP(double* v1, const double* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::min(v1[i], v2[i]);
      }

      FORCE_INLINE void MaxVectorConst(double* v1, const double* v2, size_t size_in, double c)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::max(c, v2[i]);
      }

      FORCE_INLINE void MinVectorConst(double* v1, const double* v2, size_t size_in, double c)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::min(c, v2[i]);
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function subtracts two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void SubVectorsFloat(const float* v1, const float* v2, size_t size_in, float* res)
      {
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] - v2[i];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function subtracts two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void SubVectorsUFloat(const float* v1, const float* v2, size_t size_in, float* res)
      {
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] - v2[i];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates element wise product of two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultVectorsFloat(const float* v1, const float* v2, size_t size_in, float* res)
      {
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] * v2[i];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultByConstAndAddFloat(const float* v, float c, size_t size_in, float* res)
      {
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultByConstAndAddFD(const float* v, double c, size_t size_in, double* res)
      {
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultByConstAndAddDF(const double* v, double c, size_t size_in, float* res)
      {
         for (size_t l=0; l<size_in; ++l)
            res[l] += (float)(c*v[l]);
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and places result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MultByConstFloat(const float* v, float c, size_t size_in, float* res)
      {
         for (size_t l=0; l<size_in; ++l)
            res[l] = c*v[l];
      }

      FORCE_INLINE void MultByConstIPFloat(float* v, float c, size_t size_in)
      {
         for (size_t l=0; l<size_in; ++l)
            v[l] *= c;
      }

      FORCE_INLINE void AddConstIPFloat(float* v, float c, size_t size_in)
      {
         for (size_t l=0; l<size_in; ++l)
            v[l] += c;
      }

      FORCE_INLINE void AddVectorSqrFloat(float* v1, const float* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] * v2[i];
      }

      FORCE_INLINE void MaxVectorsIPFloat(float* v1, const float* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::max(v1[i], v2[i]);
      }

      FORCE_INLINE void MinVectorsIPFloat(float* v1, const float* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::min(v1[i], v2[i]);
      }

      FORCE_INLINE void MaxVectorConstFloat(float* v1, const float* v2, size_t size_in, float c)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::max(c, v2[i]);
      }

      FORCE_INLINE void MinVectorConstFloat(float* v1, const float* v2, size_t size_in, float c)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::min(c, v2[i]);
      }

      FORCE_INLINE void AddVectorSqrFD(double* v1, const float* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] * v2[i];
      }

      FORCE_INLINE void AddVectorSqrDF(float* v1, const double* v2, size_t size_in)
      {
         for (size_t i=0; i<size_in; ++i)
            v1[i] += float(v2[i] * v2[i]);
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies rectangular matrix by a vector
      // Vector dimension should be divisible by 8 (!!!)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE void MulMatrixByVectorFloat8(const float* m, const float* v, size_t n_rows, size_t n_cols, double* res)
      {
         double sum;
         for (size_t i=0; i<n_rows; ++i)
         {
            sum = 0;
            for (size_t j=0; j<n_cols; ++j)
               sum += m[i*n_cols+j]*v[j];
            res[i] = sum;
         }
      }

      FORCE_INLINE double SumVectorElements(const double* v, size_t size_in)
      {
         double val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v[i];
         return val;
      }

      FORCE_INLINE double MinVectorElement(const double* v, size_t size_in, double m = DBL_MAX)
      {
         double val = m;
         for (size_t i=0; i<size_in; ++i)
            if (val > v[i]) val = v[i];
         return val;
      }

      FORCE_INLINE double MaxVectorElement(const double* v, size_t size_in, double m = -DBL_MAX)
      {
         double val = m;
         for (size_t i=0; i<size_in; ++i)
            if (val < v[i]) val = v[i];
         return val;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a scalar product of the two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE double DotProduct(const double* v1, const double* v2, size_t size_in)
      {
         double val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return val;
      }

      FORCE_INLINE double DotProductU(const double* v1, const double* v2, size_t size_in)
      {
         double val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return val;
      }

      float SumVectorElementsFloat(const float* v, size_t size_in)
      {
         float val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v[i];
         return val;
      }

      float MinVectorElementFloat(const float* v, size_t size_in, float m = FLT_MAX)
      {
         float val = m;
         for (size_t i=0; i<size_in; ++i)
            if (val > v[i]) val = v[i];
         return val;
      }

      FORCE_INLINE float MaxVectorElementFloat(const float* v, size_t size_in, float m = -FLT_MAX)
      {
         float val = m;
         for (size_t i=0; i<size_in; ++i)
            if (val < v[i]) val = v[i];
         return val;
      }

      FORCE_INLINE float DotProductFloat(const float* v1, const float* v2, size_t size_in)
      {
         double val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return (float)val;
      }

      FORCE_INLINE float DotProductUFloat(const float* v1, const float* v2, size_t size_in)
      {
         double val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return (float)val;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a scalar product of the two float vectors
      // Result is converted to double
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      double DotProductFD(const float* v1, const float* v2, size_t size_in)
      {
         float val = 0.f;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return double(val);
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Pow/Log/Exp/Sigmoid
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      FORCE_INLINE double Pow4(double x, double lg2_of_base) {return pow(2., x*lg2_of_base); }
      FORCE_INLINE void Pow4ArrA(const double* p_x, size_t size, double* p_res, double lg2_of_base)
      {
         for (size_t i=0; i<size; ++i)
            p_res[i] = pow(2., p_x[i]*lg2_of_base);
      }

      FORCE_INLINE double Log6(double x) { return log(x); }
      FORCE_INLINE void Log6ArrA(const double* p_x, size_t size, double* p_res)
      {
         for (size_t i=0; i<size; ++i)
            p_res[i] = log(p_x[i]);
      }

      FORCE_INLINE double Sigmoid4(double x) 
      {
         if (x<=-709) return 0;
         if (x>=709) return 1;
         
         return 1./(1.+exp(-x));
      }
      FORCE_INLINE double Sigmoid4Neg(double x)
      {
         if (x>=709) return 0;
         if (x<=-709) return 1;
         
         return 1./(1.+exp(x));
      }

      FORCE_INLINE void Sigmoid4ArrA(const double* p_x, size_t size, double* p_res)
      {
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-709) p_res[i] = 0;
            else if (p_x[i]>=709) p_res[i] = 1;
            else p_res[i] = 1./(1.+exp(-p_x[i]));
         }
      }

      FORCE_INLINE void Sigmoid4NegArrA(const double* p_x, size_t size, double* p_res)
      {
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-709) p_res[i] = 1;
            else if (p_x[i]>=709) p_res[i] = 0;
            else p_res[i] = 1./(1.+exp(p_x[i]));
         }
      }

      FORCE_INLINE double LogSigmoidL6E4(double x)
      {
         if (x<=-709) return x;
         if (x>=709) return 0;
         return -log(1.+exp(-x));
      }
      FORCE_INLINE double LogSigmoidL6E4Neg(double x)
      {
         if (x<=-709) return 0;
         if (x>=709) return -x;
         return -log(1.+exp(x));
      }
      FORCE_INLINE float LogSigmoidL6E4Neg(float x)
      {
         if (x<=-709) return 0;
         if (x>=709) return -x;
         return -log(1.f+exp(x));
      }

      FORCE_INLINE void LogSigmoidL6E4ArrA(const double* p_x, size_t size, double* p_res)
      {
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-709) p_res[i] = p_x[i];
            else if (p_x[i]>=709) p_res[i] = 0;
            else p_res[i] = -log(1.+exp(-p_x[i]));
         }
      }
      FORCE_INLINE void LogSigmoidL6E4NegArrA(const double* p_x, size_t size, double* p_res)
      {
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-709) p_res[i] = 0;
            else if (p_x[i]>=709) p_res[i] = -p_x[i];
            else p_res[i] = -log(1.+exp(p_x[i]));
         }
      }

      FORCE_INLINE float Pow4Float(float x, float lg2_of_base) {return pow(2.f, x*lg2_of_base); }
      FORCE_INLINE void Pow4ArrAFloat(const float* p_x, size_t size, float* p_res, float lg2_of_base)
      {
         for (size_t i=0; i<size; ++i)
            p_res[i] = pow(2.f, p_x[i]*lg2_of_base);
      }

      FORCE_INLINE float Log6Float(float x) { return log(x); }
      FORCE_INLINE void Log6ArrAFloat(const float* p_x, size_t size, float* p_res)
      {
         for (size_t i=0; i<size; ++i)
            p_res[i] = log(p_x[i]);
      }

      FORCE_INLINE float LogSigmoidL6E4Float(float x)
      {
         if (x<=-88) return x;
         if (x>=88) return 0;
         return -log(1.f+exp(-x));
      }
      FORCE_INLINE float LogSigmoidL6E4NegFloat(float x)
      {
         if (x<=-88) return 0;
         if (x>=88) return -x;
         return -log(1.f+exp(x));
      }

      FORCE_INLINE void LogSigmoidL6E4ArrAFloat(const float* p_x, size_t size, float* p_res)
      {
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-88) p_res[i] = p_x[i];
            else if (p_x[i]>=88) p_res[i] = 0;
            else p_res[i] = -log(1.f+exp(-p_x[i]));
         }
      }
      FORCE_INLINE void LogSigmoidL6E4NegArrAFloat(const float* p_x, size_t size, float* p_res)
      {
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-88) p_res[i] = 0;
            else if (p_x[i]>=88) p_res[i] = -p_x[i];
            else p_res[i] = -log(1.f+exp(p_x[i]));
         }
      }

#  else

      extern "C" double CityBlockDist(const double* v1, const double* v2, size_t size_in);
      /*{
         double res = 0;
         for (size_t i=0; i<size_in; ++i)
            res += fabs(v1[i]-v2[i]);
         return res;
      }*/

      extern "C" double EuclideanDist(const double* v1, const double* v2, size_t size_in);
      /*{
         double ret_val = 0;
         for (size_t i=0; i<size_in; ++i)
            ret_val += ns_base::sqr(v1[i]-v2[i]);
         return ret_val;
      }*/

      extern "C" float EuclideanDistFloat(const float* v1, const float* v2, size_t size_in);
      /*{
         float ret_val = 0;
         for (size_t i=0; i<size_in; ++i)
            ret_val += ns_base::sqr(v1[i]-v2[i]);
         return ret_val;
      }*/

      extern "C" void MulSquareMatrixByVector(const double* m, const double* v, size_t size_in, double* res);      
      /*{
         double val;
         for (size_t i=0; i<size_in; ++i)
         {
            val = 0;
            for (size_t j=0; j<size_in; ++j)
               val += m[i*size_in+j]*v[j];
            res[i] = val;
         }
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies rectangular matrix by a vector
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MulMatrixByVector(const double* m, const double* v, size_t n_rows, size_t n_cols, double* res);
      /*{
         double sum;
         for (size_t i=0; i<n_rows; ++i)
         {
            sum = 0;
            for (size_t j=0; j<n_cols; ++j)
               sum += m[i*n_cols+j]*v[j];
            res[i] = sum;
         }
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies rectangular matrix by a vector
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MulMatrixByVectorDF(const double* m, const double* v, size_t n_rows, size_t n_cols, float* res);
      /*{
         double sum;
         for (size_t i=0; i<n_rows; ++i)
         {
            sum = 0;
            for (size_t j=0; j<n_cols; ++j)
               sum += m[i*n_cols+j]*v[j];
            res[i] = sum;
         }
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies rectangular matrix by a vector
      // Vector dimension should be divisible by 8
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MulMatrixByVectorFloat(const float* m, const float* v, size_t n_rows, size_t n_cols, float* res);
      /*{
         float sum;
         for (size_t i=0; i<n_rows; ++i)
         {
            sum = 0;
            for (size_t j=0; j<n_cols; ++j)
               sum += m[i*n_cols+j]*v[j];
            res[i] = sum;
         }
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies rectangular matrix by a vector
      // Vector dimension should be divisible by 8
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MulMatrixByVectorFloat8(const float* m, const float* v, size_t n_rows, size_t n_cols, double* res);
      /*{
         double sum;
         for (size_t i=0; i<n_rows; ++i)
         {
            sum = 0;
            for (size_t j=0; j<n_cols; ++j)
               sum += m[i*n_cols+j]*v[j];
            res[i] = sum;
         }
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function is multiplying upper-triangular matrix by vector. Matrix and vector should be padded with 0 to the multiple of 4
      // Size_in is PADDED SIZE - not the original one
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MulTriangularMatrixByVector4(const double* m, const double* v, size_t size_in, double* res);
      /*{
         for (size_t i=0; i<size_in; ++i)
         {
            double tmp = 0;
            for (size_t j=i; j<size_in; ++j)
               tmp += v[j]*m[j];
            m += size_in;
            res[i] = tmp;
         }
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function subtracts two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void SubVectors(const double* v1, const double* v2, size_t size_in, double* res);
      /*{
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] - v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function subtracts two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void SubVectorsU(const double* v1, const double* v2, size_t size_in, double* res);
      /*{
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] - v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates element wise product of two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultVectors(const double* v1, const double* v2, size_t size_in, double* res);
      /*{
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] * v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorsIP(double* v1, const double* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second and the third vectors to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorsIP3(double* v1, const double* v2, const double* v3, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] + v3[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second, the third and the fourth vectors to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorsIP4(double* v1, const double* v2, const double* v3, const double* v4, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] + v3[i] + v4[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConstAndAdd(const double* v, double c, size_t size_in, double* res);
      /*{
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }*/
      
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and places result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConst(const double* v, double c, size_t size_in, double* res);
      /*{
         for (size_t l=0; l<size_in; ++l)
            res[l] = c*v[l];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant in place
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConstIP(double* v, double c, size_t size_in);
      /*{
         for (size_t l=0; l<size_in; ++l)
            v[l] *= c;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds constant to vectors in place
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddConstIP(double* v, double c, size_t size_in);
      /*{
         for (size_t l=0; l<size_in; ++l)
            v[l] += c;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorSqr(double* v1, const double* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] * v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MaxVectorsIP(double* v1, const double* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::max(v1[i], v2[i]);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MinVectorsIP(double* v1, const double* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::min(v1[i], v2[i]);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MaxVectorConst(double* v1, const double* v2, size_t size_in, double c);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::max(c, v2[i]);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MinVectorConst(double* v1, const double* v2, size_t size_in, double c);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::min(c, v2[i]);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function subtracts two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void SubVectorsFloat(const float* v1, const float* v2, size_t size_in, float* res);
      /*{
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] - v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function subtracts two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void SubVectorsUFloat(const float* v1, const float* v2, size_t size_in, float* res);
      /*{
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] - v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates element wise product of two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultVectorsFloat(const float* v1, const float* v2, size_t size_in, float* res);
      /*{
         for (size_t i=0; i<size_in; ++i)
            res[i] = v1[i] * v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorsIPFloat(float* v1, const float* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second and the third vectors to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorsIP3Float(float* v1, const float* v2, const float* v3, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] + v3[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second, the third and the fourth vectors to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorsIP4Float(float* v1, const float* v2, const float* v3, const float* v4, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] + v3[i] + v4[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConstAndAddFloat(const float* v, float c, size_t size_in, float* res);
      /*{
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConstAndAddFD(const float* v, double c, size_t size_in, double* res);
      /*{
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConstAndAddDF(const double* v, double c, size_t size_in, float* res);
      /*{
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and places result to the other vector
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConstFloat(const float* v, float c, size_t size_in, float* res);
      /*{
         for (size_t l=0; l<size_in; ++l)
            res[l] = c*v[l];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant in place
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConstIPFloat(float* v, float c, size_t size_in);
      /*{
         for (size_t l=0; l<size_in; ++l)
            v[l] *= c;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds constant to vectors in place
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddConstIPFloat(float* v, float c, size_t size_in);
      /*{
         for (size_t l=0; l<size_in; ++l)
            v[l] += c;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorSqrFloat(float* v1, const float* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] * v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MaxVectorsIPFloat(float* v1, const float* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::max(v1[i], v2[i]);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MinVectorsIPFloat(float* v1, const float* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::min(v1[i], v2[i]);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MaxVectorConstFloat(float* v1, const float* v2, size_t size_in, float c);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::max(c, v2[i]);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MinVectorConstFloat(float* v1, const float* v2, size_t size_in, float c);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] = std::min(c, v2[i]);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorSqrFD(double* v1, const float* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] * v2[i];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function adds the second vector squared element wise to the first one
      // Everything should be aligned to at least 8 (16 is prefered)
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void AddVectorSqrDF(float* v1, const double* v2, size_t size_in);
      /*{
         for (size_t i=0; i<size_in; ++i)
            v1[i] += v2[i] * v2[i];
      }*/

	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function multiplies vectors by a constant and adds result to the other vector
      // size must be devisible by 4 and everything should be aligned to 16 bytes!!!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" void MultByConstAndAddFloat4(const float* v, float c, size_t size_in, float* res);
      /*{
         for (size_t l=0; l<size_in; ++l)
            res[l] += c*v[l];
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a scalar product of the two vectors
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" double DotProduct(const double* v1, const double* v2, size_t size_in);
      /*{
         double val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a scalar product of the two vectors, unaligned version
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" double DotProductU(const double* v1, const double* v2, size_t size_in);
      /*{
         double val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a sum of vector elements
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" double SumVectorElements(const double* v, size_t size_in);
      /*{
         double val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a minimal vector element
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" double MinVectorElement(const double* v, size_t size_in, double m = DBL_MAX);
      /*{
         double val = m;
         for (size_t i=0; i<size_in; ++i)
            if (val > v[i]) val = v[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a maximal vector element
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" double MaxVectorElement(const double* v, size_t size_in, double m = -DBL_MAX);
      /*{
         double val = m;
         for (size_t i=0; i<size_in; ++i)
            if (val < v[i]) val = v[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a scalar product of the two float vectors
      // Result is float
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" float DotProductFloat(const float* v1, const float* v2, size_t size_in);
      /*{
         float val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a scalar product of the two float vectors, unaligned version
      // Result is float
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" float DotProductUFloat(const float* v1, const float* v2, size_t size_in);
      /*{
         float val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a sum of vector elements
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" float SumVectorElementsFloat(const float* v, size_t size_in);
      /*{
         float val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a minimal vector element
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" float MinVectorElementFloat(const float* v, size_t size_in, float m = FLT_MAX);
      /*{
         float val = m;
         for (size_t i=0; i<size_in; ++i)
            if (val > v[i]) val = v[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a maximal vector element
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" float MaxVectorElementFloat(const float* v, size_t size_in, float m = -FLT_MAX);
      /*{
         float val = m;
         for (size_t i=0; i<size_in; ++i)
            if (val < v[i]) val = v[i];
         return val;
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Function calculates a scalar product of the two float vectors
      // Result is converted to double
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" double DotProductFD(const float* v1, const float* v2, size_t size_in);
      /*{
         float val = 0;
         for (size_t i=0; i<size_in; ++i)
            val += v1[i]*v2[i];
         return double(val);
      }*/

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Pow/Log/Exp/Sigmoid
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      extern "C" double Pow4(double x, double lg2_of_base);
      //{return pow(2., x*lg2_of_base); }
      extern "C" void Pow4ArrA(const double* p_x, size_t size, double* p_res, double lg2_of_base);
      /*{
         for (size_t i=0; i<size; ++i)
            p_res[i] = pow(2., p_x[i]*lg2_of_base);
      }*/

      extern "C" double Log6(double x);
      //{ return log(x); }
      extern "C" void Log6ArrA(const double* p_x, size_t size, double* p_res);
      /*{
         for (size_t i=0; i<size; ++i)
            p_res[i] = log(p_x[i]);
      }*/

      extern "C" double Sigmoid4(double x);
      /*{
         if (x<=-709) return 0;
         if (x>=709) return 1;
         
         return 1./(1.+exp(-x));
      }*/
      extern "C" double Sigmoid4Neg(double x);
      /*{
         if (x>=709) return 0;
         if (x<=-709) return 1;
         
         return 1./(1.+exp(x));
      }*/

      extern "C" void Sigmoid4ArrA(const double* p_x, size_t size, double* p_res);
      /*{
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-709) p_res[i] = 0;
            else if (p_x[i]>=709) p_res[i] = 1;
            else p_res[i] = 1./(1.+exp(-x));
         }
      }*/

      extern "C" void Sigmoid4NegArrA(const double* p_x, size_t size, double* p_res);
      /*{
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-709) p_res[i] = 1;
            else if (p_x[i]>=709) p_res[i] = 0;
            else p_res[i] = 1./(1.+exp(x));
         }
      }*/

      extern "C" double LogSigmoidL6E4(double x);
      /*{
         if (x<=-709) return x;
         if (x>=709) return 0;
         return -log(1.+exp(-x));
      }*/
      extern "C" double LogSigmoidL6E4Neg(double x);
      /*{
         if (x<=-709) return 0;
         if (x>=709) return -x;
         return -log(1.+exp(x));
      }*/

      extern "C" void LogSigmoidL6E4ArrA(const double* p_x, size_t size, double* p_res);
      /*{
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-709) p_res[i] = p_x[i];
            else if (p_x[i]>=709) p_res[i] = 0;
            else p_res[i] = -log(1.+exp(-x));
         }
      }*/
      extern "C" void LogSigmoidL6E4NegArrA(const double* p_x, size_t size, double* p_res);
      /*{
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-709) p_res[i] = 0;
            else if (p_x[i]>=709) p_res[i] = -p_x[i];
            else p_res[i] = -log(1.+exp(x));
         }
      }*/

      extern "C" float Pow4Float(float x, float lg2_of_base);
      //{return pow(2., x*lg2_of_base); }
      extern "C" void Pow4ArrAFloat(const float* p_x, size_t size, float* p_res, float lg2_of_base);
      /*{
         for (size_t i=0; i<size; ++i)
            p_res[i] = pow(2., p_x[i]*lg2_of_base);
      }*/

      extern "C" float Log6Float(float x);
      //{ return log(x); }
      extern "C" void Log6ArrAFloat(const float* p_x, size_t size, float* p_res);
      /*{
         for (size_t i=0; i<size; ++i)
            p_res[i] = log(p_x[i]);
      }*/

      extern "C" float LogSigmoidL6E4Float(float x);
      /*{
         if (x<=-88) return x;
         if (x>=88) return 0;
         return -log(1.f+exp(-x));
      }*/
      extern "C" float LogSigmoidL6E4NegFloat(float x);
      /*{
         if (x<=-88) return 0;
         if (x>=88) return -x;
         return -log(1.f+exp(x));
      }*/

      extern "C" void LogSigmoidL6E4ArrAFloat(const float* p_x, size_t size, float* p_res);
      /*{
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-88) p_res[i] = p_x[i];
            else if (p_x[i]>=88) p_res[i] = 0;
            else p_res[i] = -log(1.f+exp(-x));
         }
      }*/
      extern "C" void LogSigmoidL6E4NegArrAFloat(const float* p_x, size_t size, float* p_res);
      /*{
         for (size_t i=0; i<size; ++i)
         {
            if (p_x[i]<=-88) p_res[i] = 0;
            else if (p_x[i]>=88) p_res[i] = -p_x[i];
            else p_res[i] = -log(1.f+exp(x));
         }
      }*/
#  endif

      inline double Exp4(double x) { return Pow4(x, 1.44269504088896340735992468100); }
      inline void Exp4ArrA(const double* p_x, size_t size, double* p_res) { Pow4ArrA(p_x, size, p_res, 1.44269504088896340735992468100); }
      inline float Exp4Float(float x) { return Pow4Float(x, 1.442695040888963f); }
      inline void Exp4ArrAFloat(const float* p_x, size_t size, float* p_res) { Pow4ArrAFloat(p_x, size, p_res, 1.442695040888963f); }
   }

   namespace ns_sse
   {
      FORCE_INLINE void MulMatrixByVectorA(const double* m, const double* v, size_t n_rows, size_t n_cols, double* res) {ns_base::MulMatrixByVector(m, v, n_rows, n_cols, res);}
      FORCE_INLINE void MulMatrixByVectorA(const float* m, const float* v, size_t n_rows, size_t n_cols, float* res) {ns_base::MulMatrixByVectorFloat(m, v, n_rows, n_cols, res);}
      FORCE_INLINE void MulMatrixByVectorA(const float* m, const float* v, size_t n_rows, size_t n_cols, double* res) {ns_base::MulMatrixByVectorFloat8(m, v, n_rows, n_cols, res);}
      FORCE_INLINE void MulMatrixByVectorA(const double* m, const double* v, size_t n_rows, size_t n_cols, float* res) {ns_base::MulMatrixByVectorDF(m, v, n_rows, n_cols, res);}

      FORCE_INLINE double EuclideanDistA(const double* v1, const double* v2, size_t size_in) {return ns_base::EuclideanDist(v1, v2, size_in);}
      FORCE_INLINE float EuclideanDistA(const float* v1, const float* v2, size_t size_in) {return ns_base::EuclideanDistFloat(v1, v2, size_in);}

      FORCE_INLINE double DotProductA(const double* v1, const double* v2, size_t size_in) {return ns_base::DotProduct(v1, v2, size_in);}
      FORCE_INLINE float DotProductA(const float* v1, const float* v2, size_t size_in) {return ns_base::DotProductFloat(v1, v2, size_in);}
      FORCE_INLINE double DotProductU(const double* v1, const double* v2, size_t size_in) {return ns_base::DotProductU(v1, v2, size_in);}
      FORCE_INLINE float DotProductU(const float* v1, const float* v2, size_t size_in) {return ns_base::DotProductUFloat(v1, v2, size_in);}

      FORCE_INLINE double SumVectorElementsA(const double* v, size_t size_in) {return ns_base::SumVectorElements(v, size_in);}
      FORCE_INLINE float SumVectorElementsA(const float* v, size_t size_in) {return ns_base::SumVectorElementsFloat(v, size_in);}
      FORCE_INLINE double SumVectorElementsU(const double* v, size_t size_in) {return ns_base::SumVectorElements(v, size_in);}
      FORCE_INLINE float SumVectorElementsU(const float* v, size_t size_in) {return ns_base::SumVectorElementsFloat(v, size_in);}

      FORCE_INLINE double MinVectorElementA(const double* v, size_t size_in, double m = DBL_MAX) {return ns_base::MinVectorElement(v, size_in, m);}
      FORCE_INLINE float MinVectorElementA(const float* v, size_t size_in, float m = FLT_MAX) {return ns_base::MinVectorElementFloat(v, size_in, m);}
      FORCE_INLINE double MinVectorElementU(const double* v, size_t size_in, double m = DBL_MAX) {return ns_base::MinVectorElement(v, size_in, m);}
      FORCE_INLINE float MinVectorElementU(const float* v, size_t size_in, float m = FLT_MAX) {return ns_base::MinVectorElementFloat(v, size_in, m);}

      FORCE_INLINE double MaxVectorElementA(const double* v, size_t size_in, double m = -DBL_MAX) {return ns_base::MaxVectorElement(v, size_in, m);}
      FORCE_INLINE float MaxVectorElementA(const float* v, size_t size_in, float m = -FLT_MAX) {return ns_base::MaxVectorElementFloat(v, size_in, m);}
      FORCE_INLINE double MaxVectorElementU(const double* v, size_t size_in, double m = -DBL_MAX) {return ns_base::MaxVectorElement(v, size_in, m);}
      FORCE_INLINE float MaxVectorElementU(const float* v, size_t size_in, float m = -FLT_MAX) {return ns_base::MaxVectorElementFloat(v, size_in, m);}

      FORCE_INLINE void SubVectorsA(const double* v1, const double* v2, size_t size_in, double* res) {ns_base::SubVectors(v1, v2, size_in, res);}
      FORCE_INLINE void SubVectorsA(const float* v1, const float* v2, size_t size_in, float* res) {ns_base::SubVectorsFloat(v1, v2, size_in, res);}
      FORCE_INLINE void SubVectorsU(const double* v1, const double* v2, size_t size_in, double* res) {ns_base::SubVectorsU(v1, v2, size_in, res);}
      FORCE_INLINE void SubVectorsU(const float* v1, const float* v2, size_t size_in, float* res) {ns_base::SubVectorsUFloat(v1, v2, size_in, res);}

      FORCE_INLINE void MultVectorsA(const double* v1, const double* v2, size_t size_in, double* res) {ns_base::MultVectors(v1, v2, size_in, res);}
      FORCE_INLINE void MultVectorsA(const float* v1, const float* v2, size_t size_in, float* res) {ns_base::MultVectorsFloat(v1, v2, size_in, res);}
      FORCE_INLINE void MultVectorsU(const double* v1, const double* v2, size_t size_in, double* res) {ns_base::MultVectors(v1, v2, size_in, res);}
      FORCE_INLINE void MultVectorsU(const float* v1, const float* v2, size_t size_in, float* res) {ns_base::MultVectorsFloat(v1, v2, size_in, res);}

      FORCE_INLINE void AddVectorsIPA(double* v1, const double* v2, size_t size_in) {ns_base::AddVectorsIP(v1, v2, size_in);}
      FORCE_INLINE void AddVectorsIPA(float* v1, const float* v2, size_t size_in) {ns_base::AddVectorsIPFloat(v1, v2, size_in);}
      FORCE_INLINE void AddVectorsIPU(double* v1, const double* v2, size_t size_in) {ns_base::AddVectorsIP(v1, v2, size_in);}
      FORCE_INLINE void AddVectorsIPU(float* v1, const float* v2, size_t size_in) {ns_base::AddVectorsIPFloat(v1, v2, size_in);}

      FORCE_INLINE void AddVectorsIPA(double* v1, const double* v2, const double* v3, size_t size_in) {ns_base::AddVectorsIP3(v1, v2, v3, size_in);}
      FORCE_INLINE void AddVectorsIPA(float* v1, const float* v2, const float* v3, size_t size_in) {ns_base::AddVectorsIP3Float(v1, v2, v3, size_in);}
      FORCE_INLINE void AddVectorsIPU(double* v1, const double* v2, const double* v3, size_t size_in) {ns_base::AddVectorsIP3(v1, v2, v3, size_in);}
      FORCE_INLINE void AddVectorsIPU(float* v1, const float* v2, const float* v3, size_t size_in) {ns_base::AddVectorsIP3Float(v1, v2, v3, size_in);}

      FORCE_INLINE void AddVectorsIPA(double* v1, const double* v2, const double* v3, const double* v4, size_t size_in) {ns_base::AddVectorsIP4(v1, v2, v3, v4, size_in);}
      FORCE_INLINE void AddVectorsIPA(float* v1, const float* v2, const float* v3, const float* v4, size_t size_in) {ns_base::AddVectorsIP4Float(v1, v2, v3, v4, size_in);}
      FORCE_INLINE void AddVectorsIPU(double* v1, const double* v2, const double* v3, const double* v4, size_t size_in) {ns_base::AddVectorsIP4(v1, v2, v3, v4, size_in);}
      FORCE_INLINE void AddVectorsIPU(float* v1, const float* v2, const float* v3, const float* v4, size_t size_in) {ns_base::AddVectorsIP4Float(v1, v2, v3, v4, size_in);}

      FORCE_INLINE void AddVectorSqrA(double* v1, const double* v2, size_t size_in) {ns_base::AddVectorSqr(v1, v2, size_in);}
      FORCE_INLINE void AddVectorSqrA(float* v1, const float* v2, size_t size_in) {ns_base::AddVectorSqrFloat(v1, v2, size_in);}
      FORCE_INLINE void AddVectorSqrA(double* v1, const float* v2, size_t size_in) {ns_base::AddVectorSqrFD(v1, v2, size_in);}
      FORCE_INLINE void AddVectorSqrA(float* v1, const double* v2, size_t size_in) {ns_base::AddVectorSqrDF(v1, v2, size_in);}
      FORCE_INLINE void AddVectorSqrU(double* v1, const double* v2, size_t size_in) {ns_base::AddVectorSqr(v1, v2, size_in);}
      FORCE_INLINE void AddVectorSqrU(float* v1, const float* v2, size_t size_in) {ns_base::AddVectorSqrFloat(v1, v2, size_in);}
      FORCE_INLINE void AddVectorSqrU(double* v1, const float* v2, size_t size_in) {ns_base::AddVectorSqrFD(v1, v2, size_in);}
      FORCE_INLINE void AddVectorSqrU(float* v1, const double* v2, size_t size_in) {ns_base::AddVectorSqrDF(v1, v2, size_in);}

      FORCE_INLINE void MultByConstAndAddA(const double* v, double c, size_t size_in, double* res) {ns_base::MultByConstAndAdd(v, c, size_in, res);}
      FORCE_INLINE void MultByConstAndAddA(const float* v, float c, size_t size_in, float* res) {ns_base::MultByConstAndAddFloat(v, c, size_in, res);}
      FORCE_INLINE void MultByConstAndAddA(const float* v, double c, size_t size_in, double* res) {ns_base::MultByConstAndAddFD(v, c, size_in, res);}
      FORCE_INLINE void MultByConstAndAddA(const double* v, double c, size_t size_in, float* res) {ns_base::MultByConstAndAddDF(v, c, size_in, res);}
      FORCE_INLINE void MultByConstAndAddU(const double* v, double c, size_t size_in, double* res) {ns_base::MultByConstAndAdd(v, c, size_in, res);}
      FORCE_INLINE void MultByConstAndAddU(const float* v, float c, size_t size_in, float* res) {ns_base::MultByConstAndAddFloat(v, c, size_in, res);}
      FORCE_INLINE void MultByConstAndAddU(const float* v, double c, size_t size_in, double* res) {ns_base::MultByConstAndAddFD(v, c, size_in, res);}
      FORCE_INLINE void MultByConstAndAddU(const double* v, double c, size_t size_in, float* res) {ns_base::MultByConstAndAddDF(v, c, size_in, res);}

      FORCE_INLINE void MultByConstA(const double* v, double c, size_t size_in, double* res) {ns_base::MultByConst(v, c, size_in, res);}
      FORCE_INLINE void MultByConstA(const float* v, float c, size_t size_in, float* res) {ns_base::MultByConstFloat(v, c, size_in, res);}
      FORCE_INLINE void MultByConstU(const double* v, double c, size_t size_in, double* res) {ns_base::MultByConst(v, c, size_in, res);}
      FORCE_INLINE void MultByConstU(const float* v, float c, size_t size_in, float* res) {ns_base::MultByConstFloat(v, c, size_in, res);}

      FORCE_INLINE void MultByConstIPA(double* v, double c, size_t size_in) {ns_base::MultByConstIP(v, c, size_in);}
      FORCE_INLINE void MultByConstIPA(float* v, float c, size_t size_in) {ns_base::MultByConstIPFloat(v, c, size_in);}
      FORCE_INLINE void MultByConstIPU(double* v, double c, size_t size_in) {ns_base::MultByConstIP(v, c, size_in);}
      FORCE_INLINE void MultByConstIPU(float* v, float c, size_t size_in) {ns_base::MultByConstIPFloat(v, c, size_in);}

      FORCE_INLINE void AddConstIPA(double* v, double c, size_t size_in) {ns_base::AddConstIP(v, c, size_in);}
      FORCE_INLINE void AddConstIPA(float* v, float c, size_t size_in) {ns_base::AddConstIPFloat(v, c, size_in);}
      FORCE_INLINE void AddConstIPU(double* v, double c, size_t size_in) {ns_base::AddConstIP(v, c, size_in);}
      FORCE_INLINE void AddConstIPU(float* v, float c, size_t size_in) {ns_base::AddConstIPFloat(v, c, size_in);}

      FORCE_INLINE void MaxVectorsIPA(double* v1, const double* v2, size_t size_in) {ns_base::MaxVectorsIP(v1, v2, size_in);}
      FORCE_INLINE void MaxVectorsIPA(float* v1, const float* v2, size_t size_in) {ns_base::MaxVectorsIPFloat(v1, v2, size_in);}
      FORCE_INLINE void MaxVectorsIPU(double* v1, const double* v2, size_t size_in) {ns_base::MaxVectorsIP(v1, v2, size_in);}
      FORCE_INLINE void MaxVectorsIPU(float* v1, const float* v2, size_t size_in) {ns_base::MaxVectorsIPFloat(v1, v2, size_in);}

      FORCE_INLINE void MinVectorsIPA(double* v1, const double* v2, size_t size_in) {ns_base::MinVectorsIP(v1, v2, size_in);}
      FORCE_INLINE void MinVectorsIPA(float* v1, const float* v2, size_t size_in) {ns_base::MinVectorsIPFloat(v1, v2, size_in);}
      FORCE_INLINE void MinVectorsIPU(double* v1, const double* v2, size_t size_in) {ns_base::MinVectorsIP(v1, v2, size_in);}
      FORCE_INLINE void MinVectorsIPU(float* v1, const float* v2, size_t size_in) {ns_base::MinVectorsIPFloat(v1, v2, size_in);}

      FORCE_INLINE void MaxVectorConstA(double* v1, const double* v2, size_t size_in, double c) {ns_base::MaxVectorConst(v1, v2, size_in, c);}
      FORCE_INLINE void MaxVectorConstA(float* v1, const float* v2, size_t size_in, float c) {ns_base::MaxVectorConstFloat(v1, v2, size_in, c);}
      FORCE_INLINE void MaxVectorConstU(double* v1, const double* v2, size_t size_in, double c) {ns_base::MaxVectorConst(v1, v2, size_in, c);}
      FORCE_INLINE void MaxVectorConstU(float* v1, const float* v2, size_t size_in, float c) {ns_base::MaxVectorConstFloat(v1, v2, size_in, c);}

      FORCE_INLINE void MinVectorConstA(double* v1, const double* v2, size_t size_in, double c) {ns_base::MinVectorConst(v1, v2, size_in, c);}
      FORCE_INLINE void MinVectorConstA(float* v1, const float* v2, size_t size_in, float c) {ns_base::MinVectorConstFloat(v1, v2, size_in, c);}
      FORCE_INLINE void MinVectorConstU(double* v1, const double* v2, size_t size_in, double c) {ns_base::MinVectorConst(v1, v2, size_in, c);}
      FORCE_INLINE void MinVectorConstU(float* v1, const float* v2, size_t size_in, float c) {ns_base::MinVectorConstFloat(v1, v2, size_in, c);}

      FORCE_INLINE double Pow4(double x, double lg2_of_base) {return ns_base::Pow4(x, lg2_of_base);}
      FORCE_INLINE float Pow4(float x, float lg2_of_base) {return ns_base::Pow4Float(x, lg2_of_base);}
         
      FORCE_INLINE void Pow4ArrA(const double* p_x, size_t size, double* p_res, double lg2_of_base) {return ns_base::Pow4ArrA(p_x, size, p_res, lg2_of_base);}
      FORCE_INLINE void Pow4ArrA(const float* p_x, size_t size, float* p_res, float lg2_of_base) {return ns_base::Pow4ArrAFloat(p_x, size, p_res, lg2_of_base);}

      FORCE_INLINE double Exp4(double x) {return ns_base::Exp4(x);}
      FORCE_INLINE float Exp4(float x) {return ns_base::Exp4Float(x);}
         
      FORCE_INLINE void Exp4ArrA(const double* p_x, size_t size, double* p_res) {ns_base::Exp4ArrA(p_x, size, p_res);}
      FORCE_INLINE void Exp4ArrA(const float* p_x, size_t size, float* p_res) {ns_base::Exp4ArrAFloat(p_x, size, p_res);}

      FORCE_INLINE double Log6(double x) {return ns_base::Log6(x);}
      FORCE_INLINE float Log6(float x) {return ns_base::Log6Float(x);}
         
      FORCE_INLINE void Log6ArrA(const double* p_x, size_t size, double* p_res) {ns_base::Log6ArrA(p_x, size, p_res);}
      FORCE_INLINE void Log6ArrA(const float* p_x, size_t size, float* p_res) {ns_base::Log6ArrAFloat(p_x, size, p_res);}

      FORCE_INLINE double LogSigmoidL6E4(double x) {return ns_base::LogSigmoidL6E4(x);}
      FORCE_INLINE float LogSigmoidL6E4(float x) {return ns_base::LogSigmoidL6E4Float(x);}

      FORCE_INLINE double LogSigmoidL6E4Neg(double x) {return ns_base::LogSigmoidL6E4Neg(x);}
      FORCE_INLINE float LogSigmoidL6E4Neg(float x) {return ns_base::LogSigmoidL6E4NegFloat(x);}

      FORCE_INLINE void LogSigmoidL6E4ArrA(const double* p_x, size_t size, double* p_res) {ns_base::LogSigmoidL6E4ArrA(p_x, size, p_res);}
      FORCE_INLINE void LogSigmoidL6E4ArrA(const float* p_x, size_t size, float* p_res) {ns_base::LogSigmoidL6E4ArrAFloat(p_x, size, p_res);}

      FORCE_INLINE void LogSigmoidL6E4NegArrA(const double* p_x, size_t size, double* p_res) {ns_base::LogSigmoidL6E4NegArrA(p_x, size, p_res);}
      FORCE_INLINE void LogSigmoidL6E4NegArrA(const float* p_x, size_t size, float* p_res) {ns_base::LogSigmoidL6E4NegArrAFloat(p_x, size, p_res);}
   }
}

#pragma warning(pop)

#endif