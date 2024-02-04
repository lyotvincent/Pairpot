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

const ScDendro = ({ theme, title, height, width, margin }) => {
  const chartRef = useRef(null) // get current DOM container

  useEffect(() => {
    var myChart = echarts.init(chartRef.current, theme) //init the echart container
    var _props = require('../../assets/data/Mouse-Brain-props.json')

    var cell_props = _props['pyramidal CA1']

    var line_data = []
    var lines = cell_props.lines
    for (let line of lines) {
      line_data.push({
        coords: line,
      })
    }
    myChart.setOption({
      xAxis: { type: 'value', min: 0, max: 500 },

      yAxis: { type: 'value', min: 0, max: 500, inverse: true },

      series: [
        {
          type: 'lines',
          coordinateSystem: 'cartesian2d',
          polyline: true,
          data: line_data,
          lineStyle: { color: '#00ff00', width: 5 },
          drawPoint: true,

          silent: false,
        },
      ],
    })
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

ScDendro.defaultProps = {
  theme: 'dark',
  title: 'ScDendro',
  height: '35rem',
  width: '35rem',
  margin: '2rem',
}

ScDendro.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default ScDendro
