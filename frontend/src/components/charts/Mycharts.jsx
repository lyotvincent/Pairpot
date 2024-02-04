import React, { useEffect, useRef } from 'react'
import * as echarts from 'echarts'
import '../theme/dark'

const MyCharts = (props) => {
  const chartRef = useRef() //拿到标记DOM容器，文档中是使用document获取dom。

  // 监听props，props更换时重新绘制。
  useEffect(() => {
    const chart = echarts.init(chartRef.current, 'dark') //echart初始化容器
    let option = {
      //标题
      title: {
        text: 'Referer of a Website',
        subtext: 'Fake Data',
        left: 'center',
      },
      tooltip: {
        trigger: 'item',
      },
      //图例
      legend: {
        orient: 'vertical',
        left: 'left',
      },
      series: [
        {
          name: 'Access From',
          type: 'pie',

          //可以传字符串，相当于饼图，可以传数组，相当于环状图
          radius: [20, 50],

          //一般从props中解构出Data
          data: [
            { value: 1048, name: 'Search Engine' },
            { value: 735, name: 'Direct' },
            { value: 580, name: 'Email' },
            { value: 484, name: 'Union Ads' },
            { value: 300, name: 'Video Ads' },
          ],
          emphasis: {
            itemStyle: {
              shadowBlur: 10,
              shadowOffsetX: 0,
              shadowColor: 'rgba(0, 0, 0, 0.5)',
            },
          },
        },
      ],
    }

    //将option放进 chart中
    chart.setOption(option)
  }, [props])

  return (
    <div
      ref={chartRef}
      className="chart"
      //目标dom中要记得写宽高
      style={{ height: '35rem', width: '35rem', margin: '2rem' }}></div>
  )
}

export default MyCharts
