import React, { useRef, useEffect } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import { TooltipComponent, VisualMapComponent } from 'echarts/components'
import { Scatter3DChart } from 'echarts-gl/charts'
import { Grid3DComponent } from 'echarts-gl/components'
import {
  GraphicComponent,
  GridComponent,
  LegendComponent,
} from 'echarts/components'
import { LineChart } from 'echarts/charts'
import { UniversalTransition } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
echarts.use([
  GraphicComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  LineChart,
  CanvasRenderer,
  UniversalTransition,
])

export const SpatialPlot = ({ theme, title, height, width, margin }) => {
  const chartRef = useRef() // get current DOM container

  useEffect(() => {
    var myChart = echarts.init(chartRef.current, theme) //init the echart container
    var option
    option = {
      legend: {
        data: ['Altitude (km) vs Temperature (°C)'],
      },

      tooltip: {
        trigger: 'axis',
        formatter: 'Temperature : <br/>{b}km : {c}°C',
      },

      grid: {
        left: '3%',
        right: '4%',
        bottom: '3%',
        containLabel: true,
      },

      xAxis: {
        type: 'value',
        axisLabel: {
          formatter: '{value} °C',
        },
      },

      yAxis: {
        type: 'category',
        axisLine: { onZero: false },
        axisLabel: {
          formatter: '{value} km',
        },

        boundaryGap: true,
        data: ['0', '10', '20', '30', '40', '50', '60', '70', '80'],
      },

      graphic: [
        {
          type: 'group',
          rotation: Math.PI / 4,
          bounding: 'raw',
          right: 110,
          bottom: 110,
          z: 100,
          children: [
            {
              type: 'rect',
              left: 'center',
              top: 'center',
              z: 100,
              shape: {
                width: 400,
                height: 50,
              },
              style: {
                fill: 'rgba(0,0,0,0.3)',
              },
            },
            {
              type: 'text',
              left: 'center',
              top: 'center',
              z: 100,
              style: {
                fill: '#fff',
                text: 'ECHARTS LINE CHART',
                font: 'bold 26px sans-serif',
              },
            },
          ],
        },
        {
          type: 'group',
          left: '10%',
          top: 'center',
          children: [
            {
              type: 'rect',
              z: 100,
              left: 'center',
              top: 'middle',
              shape: {
                width: 240,
                height: 90,
              },
              style: {
                fill: '#fff',
                stroke: '#555',
                lineWidth: 1,
                shadowBlur: 8,
                shadowOffsetX: 3,
                shadowOffsetY: 3,
                shadowColor: 'rgba(0,0,0,0.2)',
              },
            },
            {
              type: 'text',
              z: 100,
              left: 'center',
              top: 'middle',
              style: {
                fill: '#333',
                width: 220,
                overflow: 'break',
                text: 'xAxis represents temperature in °C, yAxis represents altitude in km, An image watermark in the upper right, This text block can be placed in any place',
                font: '14px Microsoft YaHei',
              },
            },
          ],
        },
      ],

      series: [
        {
          name: '高度(km)与气温(°C)变化关系',
          type: 'line',
          smooth: true,
          data: [15, -50, -56.5, -46.5, -22.1, -2.5, -27.7, -55.7, -76.5],
        },
      ],
    }

    option && myChart.setOption(option)
  }, [theme, title])

  return (
    <div>
      <div
        ref={chartRef}
        className="chart"
        //the target DOM container needs height and width
        style={{ height: height, width: width, margin: margin }}></div>
    </div>
  )
}

SpatialPlot.defaultProps = {
  theme: 'dark',
  title: 'SpatialPlot',
  height: '35rem',
  width: '35rem',
  margin: '2rem',
}

SpatialPlot.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default SpatialPlot
