import React, { useRef, useEffect } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import { TooltipComponent, VisualMapComponent } from 'echarts/components'
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

const CellProps = ({title, height, width, margin }) => {
  const chartRef = useRef(null) // get current DOM container

  useEffect(() => {
    var myChart = echarts.init(chartRef.current) //init the echart container
    var _props = require('../../assets/data/Mouse-Brain-props.json')
    var vega_20 = [
      '#1f77b4',
      '#aec7e8',
      '#ff7f0e',
      '#ffbb78',
      '#2ca02c',
      '#98df8a',
      '#d62728',
      '#ff9896',
      '#9467bd',
      '#c5b0d5',
      '#8c564b',
      '#c49c94',
      '#e377c2',
      '#f7b6d2',
      '#7f7f7f',
      '#c7c7c7',
      '#bcbd22',
      '#dbdb8d',
      '#17becf',
      '#9edae5',
    ]

    vega_20 = [...vega_20, ...vega_20]
    var image_row = _props.image_row
    var image_col = _props.image_col
    var cell_type = Object.keys(_props.props)
    var datasets = []
    var data_source = _props.index.map((item, idx) => {
      var ret_array = [item, image_row[idx], image_col[idx]]
      for (let ct of cell_type) {
        ret_array.push(_props.props[ct][idx])
      }
      return ret_array
    })
    datasets.push({
      dimensions: ['index', 'image_row', 'image_col', ...cell_type],
      source: data_source,
    })

    var mySeries = []
    var myVisualMap = []
    for (let id in cell_type) {
      var ct = cell_type[id]
      mySeries.push({
        type: 'scatter',
        symbolSize: 5,
        name: ct,
        encode: {
          x: 'image_row',
          y: 'image_col',
          tooltip: [0, 1, 2, ct],
        },
      })

      myVisualMap.push({
        show: ct === cell_type[0] ? true : false,
        id: ct + '_visual',
        calculable: true,
        dimension: ct,
        seriesIndex: id,
        left: 'right',
        top: '30%',
        orient: 'vertical',
        precision: 2,
        calculable: true,
        min: 0,
        max: 1,
        inRange: {
          color:["#104e8b", "#ffdab9", "#8b0a50"],
        },
        text: ['Cell-type\n proportions'],
        textGap: 20,
        textStyle: {
          fontSize: 14,
        },
      })
    }
    var line_datasets = []
    for (let id in cell_type) {
      var line_data = []
      for (let line of _props.lines[cell_type[id]]) {
        line_data.push({
          coords: line,
        })
      }
      line_datasets.push(line_data)
    }

    mySeries.push({
      type: 'lines',
      coordinateSystem: 'cartesian2d',
      polyline: true,
      data: line_datasets[0],
      lineStyle: { color: 'red', width: 5, opacity: 0.6 },
      silent: false,
      z: 10,
      visible: false,
      emphasis: {
        lineStyle: {
          color: 'green',
          width: 5,
          opacity: 1,
        },
      },
    })

    myChart.setOption({
      tooltip: {},
      color: vega_20,
      title: [
        {
          text: title,
          left: '40%',
          textStyle: {
            fontSize: 24,
          },
        },
      ],
      dataset: datasets,
      legend: {
        show: true,
        bottom: '0%',
        left: '5%',
        data: [
          ...cell_type,
          // ...cell_type.map((item) => {
          //   return item + '_lines'
          // }),
        ],
        selectedMode: 'single',
        textStyle: {
          fontSize: 16,
        },
      },
      xAxis: {
        axisTick: {
          show: false,
        },
        axisLabel: {
          show: false,
        },
        axisLine: {
          show: false,
        },
        position: 'bottom',
        min: Math.min.apply(null, image_row) - 20,
        max: Math.max.apply(null, image_row) + 20,
      },
      yAxis: {
        axisTick: {
          show: false,
        },
        axisLabel: {
          show: false,
        },
        axisLine: {
          show: false,
        },
        position: 'left',
        inverse: true,
        min: Math.min.apply(null, image_col) - 20,
        max: Math.max.apply(null, image_col) + 20,
      },
      grid: {
        top: '10%',
        left: '5%',
        right: '22%',
        bottom: '14%',
        axisLine: {
          lineStyle: {
            color: '#fff',
          },
        },
        axisPointer: {
          lineStyle: {
            color: '#ffbd67',
          },
        },
      },
      series: mySeries,
      visualMap: myVisualMap,
    })

    var currentIndex = 0

    myChart.on('legendselectchanged', (params) => {
      var isSelected = params.selected[params.name]
      var index = cell_type.indexOf(params.name)
      mySeries.pop()
      mySeries.push({
        type: 'lines',
        coordinateSystem: 'cartesian2d',
        polyline: true,
        data: line_datasets[index],
        lineStyle: { color: 'red', width: 5, opacity: 0.6 },
        silent: false,
        z: 10,
        visible: false,
        emphasis: {
          lineStyle: {
            color: 'green',
            width: 5,
            opacity: 0.8,
          },
        },
      })

      var changeVsiualMap = []
      for (let id in cell_type) {
        var ct = cell_type[id]
        changeVsiualMap.push({
          id: ct + '_visual',
          show: ct === params.name ? true : false,
        })
      }
      myChart.setOption({
        series: mySeries,
        visualMap: changeVsiualMap,
      })
      // myChart.dispatchAction({
      //   type: 'downplay',
      //   seriesIndex: currentIndex,
      // })
      // currentIndex = 2 * cell_type.indexOf(params.name) + 1
      // myChart.dispatchAction({
      //   type: 'highlight',
      //   seriesIndex: currentIndex,
      // })
    })
  }, [title])

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

CellProps.defaultProps = {
  title: 'CellProps',
  height: '35rem',
  width: '35rem',
  margin: '2rem',
}

CellProps.propTypes = {
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default CellProps
