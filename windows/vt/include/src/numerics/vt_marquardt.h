//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Class for Levenberg Marquardt minimization
//
//  History:
//      2005/7/7-swinder
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_matrix.h"
#include "vt_solve_ldl.h"

namespace vt {

#define LM_MAX_NLOOP           100
#define LM_MIN_NLOOP           1
#define LM_MIN_INCREMENT       1.0e-5
#define LM_MAX_LAMBDA          1.0e7
#define LM_LAMBDA_INCR         10.0
#define LM_LAMBDA_INIT         1.0e-3

	// Levenberg-Marquardt solver class
	//
	// Given a vector equation of the form X = f(P) where f may be a non-linear
	// function of the parameter vector P, the LM solver finds the best vector P
	// given X so as to minimise the error ||X-f(P)||.
	// LM is an extension of Newton iteration and is described in Appendix 4 of
	// Hartley Zisserman, Multiple View Geometry, and
	// http://mathworld.wolfram.com/Levenberg-MarquardtMethod.html
	//
	// Usage:
	//      Create your own class and derive it from the CSolveMarquardt class
	//      you must them inplement at least the LM_Function method to evaluate
	//      your own function. You can also supply the LM_Jacobian method if you
	//      can easily calculate the first derivative matrix, otherwise the base class
	//      calculates this derivative numerically. If an error occurs in one of
	//      your overloaded methods, then call SetError() (on CErrorBase) to prevent
	//      further iterations.

	template <class T>
	class CSolveMarquardt : public CErrorBase
	{
	public:
		CSolveMarquardt() { m_fRMSError = (T) 0.0; SetDefaults(); }
		CSolveMarquardt(const CVec<T> &initial_parameters) 
		{ SetDefaults(); SetParameters(initial_parameters); Solve(); }
		virtual ~CSolveMarquardt() {}

		void SetParameters(const CVec<T> &initial_parameters) { m_vP = initial_parameters; }

		HRESULT Solve()
		{
			CSolveLDL<T> ldl;
			CVec<T> current_cost_vector, Jte, delta, new_parameters, new_cost_vector;
			CMtx<T> J, JtJ, N;

			// This will solve the problem by Levenbert Marquardt iteration
			HRESULT hr = NOERROR;
			ClearError();

			// Get initial value of lambda
			T lambda = m_fInitialLambda;

			// Now, run for several loops
			bool terminate = false;

			VT_HR_EXIT( m_vP.GetError() );

			int loop;
			for (loop=0; loop<m_iMaximumLoops; loop++)
			{
				// The method that is used here is Gauss-Newton 
				// least squares with augmentation.  

				// First, determine the initial value of the function
				VT_HR_EXIT( LM_Function(m_vP, current_cost_vector) );
				VT_HR_EXIT( current_cost_vector.GetError() );
				VT_HR_EXIT( GetError() );

				T current_error = current_cost_vector * current_cost_vector;
				m_fRMSError = current_error;

				// Compute the Jacobian
				VT_HR_EXIT( LM_Jacobian(m_vP, J) );
				VT_HR_EXIT( J.GetError() );
				VT_HR_EXIT( GetError() );

				// We need to set up the normal equations and the goal function
				// Do it quickly to avoid taking an explicit matrix copy
				int Ndim = J.Cols();

				// Compute JtJ
				JtJ.Create(Ndim, Ndim);
				VT_HR_EXIT( JtJ.GetError() );

				int i,j;
				for (i=0; i<Ndim; i++)
					for (j=i; j<Ndim; j++)
					{
						T val = 0.0;
						for(int k=0; k<J.Rows(); k++)
							val += J[k][i] * J[k][j];
						JtJ[i][j] = val;
						JtJ[j][i] = val;
					}

					// Now, the goal vector -- actually J^T * current_cost_vector
					Jte = current_cost_vector * J;
					VT_HR_EXIT( Jte.GetError() );

					// Now, try a loop of solving and changing lambda
					for (;;)
					{
						// Have we made an improvement
						bool improvement = false; // no improvement yet

						// Now augment the normal equations
						VT_HR_EXIT( N.Create(JtJ.Rows(), JtJ.Cols()) );
						VT_HR_EXIT( LM_Augment(JtJ, lambda, N) );
						VT_HR_EXIT( N.GetError() );
						VT_HR_EXIT( GetError() );

						// Solve the Normal equations: JtJ * delta = JtJ * current_cost_vector
						delta = Jte;
						VT_HR_EXIT( delta.GetError() );

						VT_HR_EXIT( ldl.Solve(N) );
						delta = ldl.EquationSolve(delta);
						VT_HR_EXIT( delta.GetError() );

						// Try the increment
						new_parameters = m_vP - delta;
						VT_HR_EXIT( new_parameters.GetError() );

						// What is the cost
						VT_HR_EXIT( LM_Function(new_parameters, new_cost_vector) );
						VT_HR_EXIT( new_cost_vector.GetError() );
						VT_HR_EXIT( GetError() );

						T new_error = new_cost_vector * new_cost_vector;

						// Is it better than the previous?
						if (new_error < current_error)
						{
							// Yes! Decrease lambda
							lambda /= m_fLambdaIncrement;

							// Indicate that we improved
							improvement = true;
						}

						else   
						{
							// No! Increase lambda
							lambda *= m_fLambdaIncrement;

							// Indicate that we improved
							improvement = false;
						}

						// Now, determine whether to terminate or not -- before we update
						// the parameters.
						terminate = LM_CheckForTermination(m_vP, new_parameters,
							current_cost_vector, new_cost_vector,
							current_error, new_error, lambda);
						VT_HR_EXIT( GetError() );

						// Update the parameters if improved
						if (improvement)
						{
							// Accept the new parameters
							m_vP = new_parameters;

							// Give the user an opportunity to tweak the parameters.
							// For instance, it is possible to change the coordinate
							// system at this point to a local coordinate system.
							VT_HR_EXIT( LM_OnIncrement(m_vP) );
							VT_HR_EXIT( GetError() );
						}

						// Exit from this loop if terminating, or we made an improvement
						if (terminate || improvement)
							break;
					}

					// Terminate if required
					if (terminate)
						break;
			}
			m_fRMSError = sqrt(m_fRMSError/(T)current_cost_vector.Size());

Exit:
			return SetError(hr);
		}

		const CVec<T> &GetSolution() const { return m_vP; }
		const T GetRMSError() const { return m_fRMSError; }

		// These probably do not need to be used -- the defaults will do for most usage
		void SetMaximumLoops(int nloops) { m_iMaximumLoops = nloops; }
		void SetMinimumLoops(int nloops) { m_iMinimumLoops = nloops; }
		void SetMinimumIncrement(T incr) { m_fMinimumIncrement = incr; }
		void SetMaximumLambda(T lam) { m_fMaximumLambda = lam; }
		void SetLambdaIncrement(T incr) { m_fLambdaIncrement = incr; }
		void SetInitialLambda(T lam) { m_fInitialLambda = lam; }

		void SetDefaults()
		{
			m_iMaximumLoops     = LM_MAX_NLOOP; 
			m_iMinimumLoops     = LM_MIN_NLOOP; 
			m_fMinimumIncrement = (T) LM_MIN_INCREMENT; 
			m_fMaximumLambda    = (T) LM_MAX_LAMBDA;
			m_fLambdaIncrement  = (T) LM_LAMBDA_INCR;
			m_fInitialLambda    = (T) LM_LAMBDA_INIT;
		}

	protected:
		// Functions to be overloaded by derived class

		// The virtual function to be minimized -- must be supplied.
		// Note that this function must set the size of the error vector correctly
		// by using vErrorRtn.Create() as the caller does not have this information.
		virtual HRESULT LM_Function(const CVec<T> &vP, CVec<T> &vErrorRtn) = 0;

		// Same as LM_Function but called in the context of the Jacobian derivative
		// computation.
		// In that context only one parameter is changed at a time and a error delta vector
		// is computed with respect to a base value.  Super classes may over-ride
		// this if computing function in this context can be optimized.
		// This function must return f(vP+d) - f(vP). vP is the current parameter
		// set and d is the change to be added to vP. The element of vP to which d
		// must be added is indicated by iParamIndex. vErrorBase is the error vector f(vP)
		// that was previously obtained by a call to LM_Function for parameters vP.
		// This function is expected to evaluate f(vP+d) and subtract from this
		// vErrorBase returning f(vP+d) - f(vP) in vErrorDelta.
		// The element of vP indexed by iParamIndex can be modified in place without affecting
		// operation (it is auto-restored later), but no other element must be changed.
		virtual HRESULT LM_JacobianFunction(CVec<T> &vP, int iParamIndex, T d,
			const CVec<T> &vErrorBase, CVec<T> &vErrorDelta)
		{
			HRESULT hr;
			vP[iParamIndex] += d;
			VT_HR_EXIT( LM_Function(vP, vErrorDelta) );
			vErrorDelta = vErrorDelta - vErrorBase;

Exit:
			return hr;
		}

		// Operation to be carried out when an increment has occurred
		//  - default behavior is to do nothing
		virtual HRESULT LM_OnIncrement(CVec<T> &vP) { UNREFERENCED_PARAMETER(vP); 
		return NOERROR; }

		// Determine whether we should stop or continue iterating
		virtual bool LM_CheckForTermination (
			const CVec<T> &old_parameters, 
			const CVec<T> &new_parameters,
			const CVec<T> &old_cost_vector,// = function(old_parameters)
			const CVec<T> &new_cost_vector,// = function(new_parameters)
			T old_error,                   // = old_cost_vector * old_cost_vector
			T new_error,                   // = new_cost_vector * new_cost_vector
			T lambda                       // current LM lambda parameter
			) const
		{
			// Compute the change in the cost vector, and the proportional error
			int i;
			T magsq = 0;
			int size = old_cost_vector.Size();
			for(i=0; i<size; i++)
			{
				T diff = old_cost_vector[i] - new_cost_vector[i];
				magsq += diff*diff;
			}

			T error_delta = sqrt(magsq / old_error);

			// If we are taking too small steps of lambda becomes too big,
			// then we are going nowhere, so stop
			if (error_delta < m_fMinimumIncrement || lambda > m_fMaximumLambda)
				return true;
			else
				return false;
		}

		// Jacobian if supplied
		// - default behavior is to calculate the Jacobian using numerical differentiation
		// Note that this function must set the size of the J matrix correctly
		// by using mJRtn.Create() as the caller does not have this information.
		virtual HRESULT LM_Jacobian(const CVec<T> &vP, CMtx<T> &J)
		{
			HRESULT hr = NOERROR;
			CVec<T> errcostvec, costvec0;

			// Constants used in computing numberical derivatives
			T DeltaFactor = (T)1.0e-6;
			T MinDelta    = (T)1.0e-6;

			// Take a copy of the input parameters
			CVec<T> tparams = vP;
			VT_HR_EXIT( tparams.GetError() );

			// Evaluate the function at the current point
			VT_HR_EXIT( LM_Function(tparams, costvec0) );

			// Now, allocate a matrix of the correct size
			J.Create(costvec0.Size(), tparams.Size());
			VT_HR_EXIT( J.GetError() );

			// Now, vary each of the parameters, one by one 
			// and fill out the Jacobian
			int i, j;
			for (j=0; j<tparams.Size(); j++)
			{
				// Save the original value of the parameter to be changed
				T param = tparams[j];

				//find incrememnt
				T delta = (T) (param * DeltaFactor);
				if (fabs(delta) < MinDelta)
					delta = MinDelta;

				// Compute the function at this value
				VT_HR_EXIT( LM_JacobianFunction (tparams, j, delta, costvec0, errcostvec) );

				VT_ASSERT( costvec0.Size() == errcostvec.Size() );

				// Use this to fill out the derivatives
				for (i=0; i<J.Rows(); i++)
					J[i][j] = errcostvec[i] / delta;

				// restore the parameter in case it was changed
				tparams[j] = param;
			}

Exit:
			return SetError(hr);
		}


		// Augment the normal equations
		// - default behavior is to form M_ii = (1 + lambda) * M_ii
		virtual HRESULT LM_Augment (const CMtx<T> &M, const T lambda, CMtx<T> &M2)
		{
			// Return matrix with diagonals multiplied by (1+lambda)
			// HRESULT hr = M2.Create(M.Rows(),M.Cols()) );

			int min_dim = (M.Rows()<=M.Cols()) ? M.Rows() : M.Cols();
			T factor = (T)(1.0 + lambda);

			M2 = M;
			for (int i=0; i<min_dim; i++)
				M2[i][i] *= factor;

			return NOERROR;
		}

		// The current parameters
		CVec<T> m_vP;
		T m_fRMSError;

		// Execution parameters
		int m_iMaximumLoops;
		int m_iMinimumLoops;
		T m_fMinimumIncrement;
		T m_fMaximumLambda;
		T m_fLambdaIncrement;
		T m_fInitialLambda;
	};

	typedef CSolveMarquardt<float> CSolveMarquardtf;
	typedef CSolveMarquardt<double> CSolveMarquardtd;

};
