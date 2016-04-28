#pragma once

#define _USE_MATH_DEFINES
#include <comdef.h>
#include <d3d11_1.h>
#include <d2d1_1.h>
#include <wincodec.h>
#include <math.h>

#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "d2d1.lib")
#pragma comment(lib, "d2d1.lib")
#pragma comment(lib, "windowscodecs.lib")

#define SMART_PTR_DEF(x) _COM_SMARTPTR_TYPEDEF(x,__uuidof(x));
#define RETURN_ON_FAIL(x) { if (FAILED(x)) return x; }

SMART_PTR_DEF(ID3D11Device)
SMART_PTR_DEF(ID3D11DeviceContext)
SMART_PTR_DEF(ID2D1DeviceContext)
SMART_PTR_DEF(ID3D11Texture2D)
SMART_PTR_DEF(IDXGISwapChain)
SMART_PTR_DEF(IDXGISurface)
SMART_PTR_DEF(ID2D1SolidColorBrush)
SMART_PTR_DEF(ID3D10Multithread)
SMART_PTR_DEF(ID2D1Factory)
SMART_PTR_DEF(ID3D11Resource)
SMART_PTR_DEF(IWICImagingFactory)
SMART_PTR_DEF(IWICBitmapEncoder)
SMART_PTR_DEF(IWICBitmapFrameEncode)
SMART_PTR_DEF(IWICStream)
SMART_PTR_DEF(IPropertyBag2)
SMART_PTR_DEF(ID3D11Query)

class CD2D1Graph
{

public:

	CD2D1Graph(FLOAT width, FLOAT height);

	~CD2D1Graph();

	HRESULT Initialize();

	HRESULT BeginDraw();

	HRESULT DrawPointPolar(FLOAT radius, FLOAT radians);

	HRESULT DrawPoint(FLOAT x, FLOAT y);

	HRESULT EndDraw();

	HRESULT Present();

	HWND GetWindowHandle() { return m_WindowHandle; }
	
	HANDLE GetInitSyncEventHandle() { return m_InitSyncEventHandle; }

	HRESULT SavePNG(const wchar_t* fileName);
private:

	static LRESULT CALLBACK WindowProc(
		_In_ HWND   hwnd,
		_In_ UINT   uMsg,
		_In_ WPARAM wParam,
		_In_ LPARAM lParam);
	
	static DWORD WINAPI ThreadProc(
		_In_ LPVOID lpParameter
		);
	
	ID2D1FactoryPtr m_pFactory = nullptr;
	ID2D1SolidColorBrushPtr m_pSolidBrush1 = nullptr;
	ID2D1SolidColorBrushPtr m_pSolidBrush2 = nullptr;
	ID3D11DevicePtr m_pD3D11Device = nullptr;
	IDXGISwapChainPtr m_pSwapChain = nullptr;
	ID3D11DeviceContextPtr m_pD3D11DeviceContext = nullptr;
	ID2D1DeviceContextPtr m_pD2D1DeviceContext = nullptr;
	HWND m_WindowHandle;
	HANDLE m_WindowThreadHandle;
	HANDLE m_InitSyncEventHandle;
	FLOAT m_Width;
	FLOAT m_Height;

};

