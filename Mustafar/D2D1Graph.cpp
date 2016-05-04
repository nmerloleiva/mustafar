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
	CoInitialize(NULL);

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

CD2D1Graph::CD2D1Graph(FLOAT width, FLOAT height, FLOAT scale)
{
	CoInitialize(NULL);
	m_Width = width;
	m_Height = height;
	m_Scale = scale;
	m_InitSyncEventHandle = CreateEvent(NULL, FALSE, FALSE, NULL);
	m_WindowThreadHandle = CreateThread(NULL, 0, ThreadProc, this, 0, NULL);
	WaitForSingleObject(m_InitSyncEventHandle, INFINITE);
}

CD2D1Graph::~CD2D1Graph()
{
	PostMessage(GetWindowHandle(), WM_CLOSE, 0, 0);
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
	
	//ShowWindow(m_WindowHandle, 1);

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

	m_pD2D1DeviceContext->Clear(D2D1::ColorF(0, 0, 0, 0)); // Transparent background
	RETURN_ON_FAIL(hr);
	
	m_pD2D1DeviceContext->SetAntialiasMode(D2D1_ANTIALIAS_MODE_ALIASED);

	D2D1_ELLIPSE ellipse;
	ellipse.point.x = m_Width / 2;
	ellipse.point.y = m_Height / 2;
	ellipse.radiusX = 20;
	ellipse.radiusY = 20;
	m_pD2D1DeviceContext->DrawEllipse(ellipse, m_pSolidBrush1, 1.0);
	m_pD2D1DeviceContext->DrawLine(D2D1::Point2F(0, ellipse.point.y), D2D1::Point2F(m_Width, ellipse.point.y), m_pSolidBrush1);
	m_pD2D1DeviceContext->DrawLine(D2D1::Point2F(ellipse.point.x, 0), D2D1::Point2F(ellipse.point.x, m_Height), m_pSolidBrush1);

	hr = m_pD2D1DeviceContext->EndDraw();
	RETURN_ON_FAIL(hr);

	hr = m_pSwapChain->Present(0, 0);
	RETURN_ON_FAIL(hr);

	return S_OK;
}

HRESULT CD2D1Graph::BeginDraw()
{
	m_pD2D1DeviceContext->BeginDraw();
	m_pD2D1DeviceContext->SetAntialiasMode(D2D1_ANTIALIAS_MODE_PER_PRIMITIVE);
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
	D2D1_ELLIPSE ellipse;
	ellipse.point.x = x / m_Scale + (m_Width / 2);
	ellipse.point.y = y / m_Scale + (m_Height / 2);
	ellipse.radiusX = 3;
	ellipse.radiusY = 3;
	m_pD2D1DeviceContext->DrawEllipse(ellipse, m_pSolidBrush2, 1.0);
	m_pD2D1DeviceContext->FillEllipse(ellipse, m_pSolidBrush2);
	return S_OK;
}

HRESULT CD2D1Graph::EndDraw()
{
	return m_pD2D1DeviceContext->EndDraw();;
}

HRESULT CD2D1Graph::Present()
{
	return m_pSwapChain->Present(0, 0);;
}

HRESULT CD2D1Graph::SavePNG(const wchar_t* fileName)
{
	HRESULT hr = S_OK;

	ID3D11Texture2DPtr backBuffer = nullptr;
	hr = m_pSwapChain->GetBuffer(0, __uuidof(backBuffer), (void**)&backBuffer);
	RETURN_ON_FAIL(hr);

	ID3D11Texture2DPtr stagingTexture = nullptr;
	D3D11_TEXTURE2D_DESC textureDesc = { 0 };
	backBuffer->GetDesc(&textureDesc);
	textureDesc.Usage = D3D11_USAGE_STAGING;
	textureDesc.BindFlags = 0;
	textureDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;

	hr = m_pD3D11Device->CreateTexture2D(&textureDesc, nullptr, &stagingTexture);
	RETURN_ON_FAIL(hr);

	m_pD3D11DeviceContext->CopyResource(stagingTexture, backBuffer);

	IWICImagingFactoryPtr pFactory = nullptr;
	hr = CoCreateInstance(
		CLSID_WICImagingFactory,
		NULL,
		CLSCTX_INPROC_SERVER,
		IID_PPV_ARGS(&pFactory));

	IWICBitmapEncoderPtr pEncoder = nullptr;
	hr = pFactory->CreateEncoder(GUID_ContainerFormatPng, nullptr, &pEncoder);
	RETURN_ON_FAIL(hr);

	IWICStreamPtr pWICStream = nullptr;
	hr = pFactory->CreateStream(&pWICStream);
	RETURN_ON_FAIL(hr);

	hr = pWICStream->InitializeFromFilename(fileName, GENERIC_WRITE);
	RETURN_ON_FAIL(hr);

	hr = pEncoder->Initialize(pWICStream, WICBitmapEncoderNoCache);
	RETURN_ON_FAIL(hr);

	IWICBitmapFrameEncodePtr pBitmapFrameEncode = nullptr;
	hr = pEncoder->CreateNewFrame(&pBitmapFrameEncode, nullptr);
	RETURN_ON_FAIL(hr);

	hr = pBitmapFrameEncode->Initialize(nullptr);
	RETURN_ON_FAIL(hr);

	hr = pBitmapFrameEncode->SetSize(m_Width, m_Height);
	RETURN_ON_FAIL(hr);

	GUID pixelFormat = GUID_WICPixelFormat32bppBGRA;
	hr = pBitmapFrameEncode->SetPixelFormat(&pixelFormat);
	RETURN_ON_FAIL(hr);

	D3D11_MAPPED_SUBRESOURCE mapped;
	hr = m_pD3D11DeviceContext->Map(stagingTexture, 0, D3D11_MAP_READ, 0, &mapped);
	RETURN_ON_FAIL(hr);

	BYTE* pBuffer = (BYTE*)mapped.pData;
	UINT bufferSize = m_Height * mapped.RowPitch;

	hr = pBitmapFrameEncode->WritePixels(m_Height, mapped.RowPitch, bufferSize, pBuffer);
	RETURN_ON_FAIL(hr);

	hr = pBitmapFrameEncode->Commit();
	RETURN_ON_FAIL(hr);

	hr = pEncoder->Commit();
	RETURN_ON_FAIL(hr);

	m_pD3D11DeviceContext->Unmap(stagingTexture, 0);

	return S_OK;
}
