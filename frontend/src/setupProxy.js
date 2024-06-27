import createProxyMiddleware from "http-proxy-middleware"

module.exports = function (app) {
  app.use(
    createProxyMiddleware('/api', {
      target: 'http://localhost:5522',
      changeOrigin: true,
      ws: true,
      pathRewrite: { '^/api': '' }  // 这句话是将你请求中的/api换成/api
    })
  )
}
