import React, { useRef, useEffect } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import {
  GraphicComponent,
  GridComponent,
  LegendComponent,
  TooltipComponent,
  VisualMapComponent,
} from 'echarts/components'
import { LinesChart } from 'echarts/charts'
import { UniversalTransition } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
echarts.use([
  GraphicComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  CanvasRenderer,
  UniversalTransition,
  VisualMapComponent,
  LinesChart,
])

const Drag = ({ theme, title, height, width, margin }) => {
  const chartRef = useRef(null) // get current DOM container

  useEffect(() => {
    var symbolSize = 20
    var data = [
      [15, 0],
      [-50, 10],
      [-30, 15],
      [-45, 20],
      [-22, 40],
    ]
    var myChart = echarts.init(chartRef.current, theme) //init the echart container

    myChart.setOption({
      xAxis: { type: 'value', min: -100, max: 80, axisLine: { onZero: false } },

      yAxis: {
        type: 'value',
        min: -30,
        max: 60,
        axisLine: { onZero: false },
      },

      series: [
        {
          type: 'line',
          id: 'a',
          smooth: true,
          symbolSize: symbolSize,
          data: data,
        },
      ],
    })

    myChart.setOption({
      graphic: echarts.util.map(data, function (dataItem, dataIndex) {
        return {
          type: 'circle',
          shape: {
            r: symbolSize / 2,
          },

          position: myChart.convertToPixel('grid', dataItem),

          invisible: true,
          draggable: true,
          z: 100,
          ondrag: echarts.util.curry(onPointDragging, dataIndex),
        }
      }),
    })

    function onPointDragging(dataIndex, dx, dy) {
      data[dataIndex] = myChart.convertFromPixel('grid', this.position)
      myChart.setOption({
        series: [
          {
            id: 'a',
            data: data,
          },
        ],
      })
    }
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

Drag.defaultProps = {
  theme: 'dark',
  title: 'Drag',
  height: '35rem',
  width: '35rem',
  margin: '2rem',
}

Drag.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default Drag
