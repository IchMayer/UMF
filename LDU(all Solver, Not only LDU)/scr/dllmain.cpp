// dllmain.cpp : Определяет точку входа для приложения DLL.
#include "stdafx.h"
#include "LDU.cpp"

Matrix::Gause::Gause<double> gause;
Matrix::LDU::LDU<double> ldu;

extern "C"
{
	_declspec(dllexport) void Gause(std::vector<std::vector<double>> &A, std::vector<double> &F, std::vector<double> &X)
	{
		gause.A = A;
		gause.f = F;
		Matrix::Gause::CountUpTrianMatrix(gause);
		Matrix::Gause::GausBack(gause);
		X = gause.x;
	}

	_declspec(dllexport) void LDU(std::vector<double> al, std::vector<double> au, std::vector<double> di,  std::vector<size_t> &ia, std::vector<double> &F, std::vector<double> &X)
	{
		ldu.n = di.size();
		ldu.m = al.size();
		ldu.x = X;

		Matrix::LDU::CountLDU(ldu, al, au, di, ia);
		Matrix::LDU::CountX(ldu, F);

		X = ldu.x;
	}

	_declspec(dllexport) size_t MSGSolver(std::vector<size_t> &iig,std::vector<size_t> &ijg,std::vector<double> &gl,std::vector<double> &gu,std::vector<double> &di, std::vector<double> &F, std::vector<double> &X, double EPS)
	{
		Matrix::MSG_LOS_BCGSTAB::Matrix<double> msg;
		size_t n, m;
		n = di.size();
		m = gl.size();
		msg.Init(n, m, iig, ijg, gl, gu, di, F, EPS);

		msg.preconditioning(3);
		msg.ConjugateGradient(3);

		X = msg.Result();
		return msg.K();
	}

	_declspec(dllexport) size_t LOS(std::vector<size_t> &iig, std::vector<size_t> &ijg, std::vector<double> &gl, std::vector<double> &gu, std::vector<double> &di, std::vector<double> &F, std::vector<double> &X, double EPS)
	{
		Matrix::MSG_LOS_BCGSTAB::Matrix<double> msg;
		size_t n, m;
		n = di.size();
		m = gl.size();
		msg.Init(n, m, iig, ijg, gl, gu, di, F, EPS);

		msg.preconditioning(3);
		msg.LocalOptimalScheme(3);

		X = msg.Result();

		return msg.K();
	}

	_declspec(dllexport) size_t BSG_STAB(std::vector<size_t> &iig, std::vector<size_t> &ijg, std::vector<double> &gl, std::vector<double> &gu, std::vector<double> &di, std::vector<double> &F, std::vector<double> &X, double EPS)
	{
		Matrix::MSG_LOS_BCGSTAB::Matrix<double> msg;
		size_t n, m;
		n = di.size();
		m = gl.size();
		msg.Init(n, m, iig, ijg, gl, gu, di, F, EPS);

		msg.preconditioning(3);
		msg.BSG_STAB(3);

		X = msg.Result();

		return msg.K();
	}
}

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
	setlocale(LC_ALL, "Russian");
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

