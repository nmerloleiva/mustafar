#include "stdafx.h"
#include "D2D1Graph.h"

LRESULT CALLBACK CD2D1Graph::WindowProc(
	_In_ HWND   hwnd,
	_In_ UINT   uMsg,
	_In_ WPARAM wParam,
	_In_ LPARAM lParam)
{

	return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

DWORD WINAPI CD2D1Graph::ThreadProc(
	_In_ LPVOID lpParameter)
{
	CD2D1Graph* graph = (CD2D1Graph*)lpParameter;
	graph->Initialize();
	SetEvent(graph->GetInitSyncEventHandle());

	HWND windowHandle = graph->GetWindowHandle();
	MSG msg;
	while (GetMessage(&msg, windowHandle, NULL, NULL) > 0){
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return 0;
}

CD2D1Graph::CD2D1Graph(FLOAT width, FLOAT height)
{
	m_Width = width;
	m_Height = height;
	m_InitSyncEventHandle = CreateEvent(NULL, FALSE, FALSE, NULL);
	m_WindowThreadHandle = CreateThread(NULL, 0, ThreadProc, this, 0, NULL);
	WaitForSingleObject(m_InitSyncEventHandle, INFINITE);
}

CD2D1Graph::~CD2D1Graph()
{
	WaitForSingleObject(m_WindowThreadHandle, INFINITE);
}

HRESULT CD2D1Graph::Initialize()
{
	HRESULT hr = S_OK;

	WNDCLASSEX wndclassex = { 0 };
	wndclassex.cbSize = sizeof(WNDCLASSEX);
	wndclassex.hInstance = GetModuleHandle(NULL);
	wndclassex.lpfnWndProc = WindowProc;
	wndclassex.lpszClassName = L"MustafarWindow";
	RegisterClassEx(&wndclassex);

	m_WindowHandle = CreateWindowEx(WS_EX_COMPOSITED, L"MustafarWindow", L"Mustafar", 0, 0, 0, m_Width, m_Height, NULL, NULL, 0, this);
	
	ShowWindow(m_WindowHandle, 1);

	D3D_FEATURE_LEVEL featureLevelRequest[] = { D3D_FEATURE_LEVEL_11_1, D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_1, D3D_FEATURE_LEVEL_10_0 };
	UINT featureLevelRequestCount = sizeof(featureLevelRequest) / sizeof(featureLevelRequest[0]);

	D3D_FEATURE_LEVEL featureLevel;
	
	DXGI_SWAP_CHAIN_DESC swapChainDesc = { 0 };
	swapChainDesc.BufferDesc.Width = m_Width;
	swapChainDesc.BufferDesc.Height = m_Height;
	swapChainDesc.BufferDesc.RefreshRate = {60, 1};
	swapChainDesc.BufferDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM,
	swapChainDesc.BufferDesc.ScanlineOrdering = DXGI_MODE_SCANLINE_ORDER_PROGRESSIVE,
	swapChainDesc.BufferDesc.Scaling = DXGI_MODE_SCALING_UNSPECIFIED;
	swapChainDesc.SampleDesc.Count = 1;
	swapChainDesc.SampleDesc.Quality = 0;
	swapChainDesc.BufferCount = 2;
	swapChainDesc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	swapChainDesc.OutputWindow = m_WindowHandle;
	swapChainDesc.Windowed = TRUE;
	swapChainDesc.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;
	swapChainDesc.Flags = 0;

	hr = D3D11CreateDeviceAndSwapChain(
		nullptr,
		D3D_DRIVER_TYPE_HARDWARE,
		NULL,
		D3D11_CREATE_DEVICE_BGRA_SUPPORT,
		featureLevelRequest,
		featureLevelRequestCount,
		D3D11_SDK_VERSION,
		&swapChainDesc,
		&m_pSwapChain,
		&m_pD3D11Device,
		&featureLevel,
		&m_pD3D11DeviceContext);
	RETURN_ON_FAIL(hr);

	hr = ((ID3D10MultithreadPtr)m_pD3D11DeviceContext)->SetMultithreadProtected(TRUE);
	RETURN_ON_FAIL(hr);

	IDXGISurfacePtr backBuffer = nullptr;
	hr = m_pSwapChain->GetBuffer(0, __uuidof(backBuffer), (void**)&backBuffer);
	RETURN_ON_FAIL(hr);

	D2D1_CREATION_PROPERTIES deviceContextProps = D2D1::CreationProperties(D2D1_THREADING_MODE_MULTI_THREADED, D2D1_DEBUG_LEVEL_NONE, D2D1_DEVICE_CONTEXT_OPTIONS_NONE);
	hr = D2D1CreateDeviceContext(backBuffer, deviceContextProps, &m_pD2D1DeviceContext);
	RETURN_ON_FAIL(hr);

	hr = D2D1CreateFactory(D2D1_FACTORY_TYPE_MULTI_THREADED, &m_pFactory);
	RETURN_ON_FAIL(hr);

	FLOAT dpiX, dpiY;
	m_pFactory->GetDesktopDpi(&dpiX, &dpiY);
	m_pD2D1DeviceContext->SetDpi(dpiX, dpiY);

	hr = m_pD2D1DeviceContext->CreateSolidColorBrush(D2D1::ColorF(0, 0, 0, 1), &m_pSolidBrush1);
	RETURN_ON_FAIL(hr);

	hr = m_pD2D1DeviceContext->CreateSolidColorBrush(D2D1::ColorF(0, 0, 0, 1), &m_pSolidBrush2);
	RETURN_ON_FAIL(hr);

	m_pD2D1DeviceContext->BeginDraw();

	m_pD2D1DeviceContext->Clear(D2D1::ColorF(1, 1, 1, 1));
	RETURN_ON_FAIL(hr);
	
	m_pD2D1DeviceContext->SetAntialiasMode(D2D1_ANTIALIAS_MODE_ALIASED);

	m_pD2D1DeviceContext->DrawLine(D2D1::Point2F(m_Width / 2 - 5, m_Height / 2), D2D1::Point2F(m_Width / 2 + 5, m_Height / 2), m_pSolidBrush1);

	m_pD2D1DeviceContext->DrawLine(D2D1::Point2F(m_Width / 2, m_Height / 2 - 5), D2D1::Point2F(m_Width / 2, m_Height / 2 + 5), m_pSolidBrush1);

	hr = m_pD2D1DeviceContext->EndDraw();
	RETURN_ON_FAIL(hr);

	hr = m_pSwapChain->Present(0, 0);
	RETURN_ON_FAIL(hr);

	return S_OK;
}

HRESULT CD2D1Graph::BeginDraw()
{
	m_pD2D1DeviceContext->BeginDraw();
	m_pD2D1DeviceContext->SetAntialiasMode(D2D1_ANTIALIAS_MODE_ALIASED);
	return S_OK;
}

HRESULT CD2D1Graph::DrawPointPolar(FLOAT radius, FLOAT radians)
{
	FLOAT x = radius * cos(radians);
	FLOAT y = radius * sin(radians);
	return DrawPoint(x, y);
}

HRESULT CD2D1Graph::DrawPoint(FLOAT x, FLOAT y)
{
	D2D1_RECT_F rect;
	rect.left = rect.right = x + (m_Width / 2);
	rect.top = rect.bottom = -y + (m_Height / 2);
	m_pD2D1DeviceContext->DrawRectangle(rect, m_pSolidBrush2, 1.0);
	return S_OK;
}

HRESULT CD2D1Graph::EndDraw()
{
	return  m_pD2D1DeviceContext->EndDraw();;
}

HRESULT CD2D1Graph::Present()
{
	return m_pSwapChain->Present(0, 0);;
}