import React, { useRef, useEffect, useState } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import { TooltipComponent, VisualMapComponent } from 'echarts/components'
import { CanvasRenderer } from 'echarts/renderers'
import { Scatter3DChart } from 'echarts-gl/charts'
import { Grid3DComponent } from 'echarts-gl/components'

echarts.use([
  TooltipComponent,
  VisualMapComponent,
  Grid3DComponent,
  Scatter3DChart,
  CanvasRenderer,
])

export const Landscape3D2 = ({
  title,
  height,
  width,
  margin,
  initialVisible = true,
}) => {
  const [visible, setVisible] = useState(initialVisible)
  const chartRef = useRef() // get current DOM container

  useEffect(() => {
    var app = {}
    var myChart = echarts.init(chartRef.current) //init the echart container
    var option
    var data
    var _data = require('../../assets/data/PDAC-3D.json')
    var cell_types = []
    var colors = []

    var datasets = []
    datasets.push({
      dimensions: ['X', 'Y', 'Z', 'color', 'cell_type'],
      source: _data.map((item) => {
        return [
          item['X'],
          item['Y'],
          item['Z'],
          item['color'],
          item['cell_type'],
        ]
      }),
    })
    var currentGroup
    var seriesArray = []
    var seriesNum = 0
    for (var i = 0; i < _data.length; i++) {
      // seperate the _data into subgroups by cell-types
      var item = _data[i]
      if (!currentGroup || currentGroup.cell_type !== item.cell_type) {
        seriesNum = seriesNum + 1
        currentGroup = { cell_type: item.cell_type }
        cell_types.push(item.cell_type)
        colors.push(item.color)
        datasets.push({
          // 这个 dataset 的 index 是 `1`。
          transform: {
            type: 'filter',
            config: { dimension: 'cell_type', value: item.cell_type },
          },
        })
        seriesArray.push({
          type: 'scatter3D',
          name: item.cell_type,
          encode: {
            x: 'X',
            y: 'Y',
            z: 'Z',
            tooltip: [0, 1, 4],
            itemName: item.cell_type,
          },
          datasetIndex: seriesNum,
        })
      }
    }

    var VisualColors = JSON.parse(JSON.stringify(colors)).reverse()
    myChart.setOption({
      title: [
        {
          text: title,
          left: '47.5%',
          textStyle: {
            fontSize: 24,
          },
        },
      ],
      color: colors,
      tooltip: {},
      xAxis3D: {},
      yAxis3D: {},
      zAxis3D: {},
      grid3D: {
        left: '-15%',
        top: '-5%',
        // axisLine: {
        //   lineStyle: {
        //     color: '#fff',
        //   },
        // },
        // axisPointer: {
        //   lineStyle: {
        //     color: '#ffbd67',
        //   },
        // },
        viewControl: {
          autoRotate: true,
          // projection: 'orthographic'
        },
      },
      legend: {
        orient: 'vertical',
        left: '70%',
        top: 'center',
        textStyle: {
          fontSize: 16,
        },
      },
      toolbox: {
        show: true,
        orient: 'horizontal',
        left: 'right',
        top: '85%',
        feature: {
          mark: { show: true },
          dataView: { show: true, readOnly: false },
          restore: { show: true },
          saveAsImage: { show: true },
        },
      },
      dataset: datasets,
      visualMap: {
        show: false,
        type: 'piecewise',
        dimension: 'cell_type',
        categories: cell_types,
        color: VisualColors,
      },
      series: seriesArray,
    })
  }, [ title, visible])

  const handleToggle = () => {
    setVisible(!visible)
  }

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

Landscape3D2.defaultProps = {
  title: 'Landscape3D2',
  height: '35rem',
  width: '70rem',
  margin: '2rem',
}

Landscape3D2.propTypes = {
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default Landscape3D2
