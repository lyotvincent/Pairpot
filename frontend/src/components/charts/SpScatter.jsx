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

const SpScatter = ({ theme, title, height, width, margin }) => {
  const chartRef = useRef(null) // get current DOM container

  useEffect(() => {
    var myChart = echarts.init(chartRef.current, theme) //init the echart container
    var _data = require('../../assets/data/Mouse-Brain-umap.json')

    var array_row = _data.map((item) => {
      return -item['array_row']
    })
    var array_col = _data.map((item) => {
      return -item['array_col']
    })

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

    var umapSeries = [
      {
        type: 'scatter',
        symbolSize: 6,
        name: 'V1_Adult_Mouse_Brain',
        dimensions: ['array_col', 'array_row', 'umap_1', 'umap_2', 'umap_3'],
        data: _data.map((item) => {
          return [
            item['array_col'],
            -item['array_row'],
            item['umap_1'],
            item['umap_2'],
            item['umap_3'],
          ]
        }),
        encode: {
          x: 'array_col',
          y: 'array_row',
          tooltip: [0, 1, 2, 3, 4],
        },
        itemStyle: {
          color: (params) => {
            return (
              'rgba(' +
              params.data[2] * 255 +
              ',' +
              params.data[3] * 255 +
              ',' +
              params.data[4] * 255 +
              ',1)'
            )
          },
        },
      },
    ]
    var legendData = ['V1_Adult_Mouse_Brain']
    var cluSeries = []
    var item = _data[0]
    for (let key of Object.keys(item)) {
      if (
        [
          'BayesSpace',
          'SEDR',
          'stlearn',
          'Giotto',
          'Seurat',
          'leiden',
          'SpaGCN',
        ].includes(key)
      ) {
        legendData.push(key)
        cluSeries.push({
          type: 'scatter',
          symbolSize: 6,
          name: key,
          dimensions: ['array_col', 'array_row', 'Cluster'],
          data: _data.map((item) => {
            return [item['array_col'], -item['array_row'], item[key]]
          }),
          encode: {
            x: 'array_col',
            y: 'array_row',
            tooltip: [0, 1, 2],
          },
          itemStyle: {
            color: (params) => {
              return vega_20[params.data[2]]
            },
          },
        })
      }
    }
    var mySeries = [...umapSeries, ...cluSeries]
    myChart.setOption({
      tooltip: {},
      title: [
        {
          text: title,
          left: 'center',
          textStyle: {
            fontSize: 24,
          },
        },
      ],
      xAxis: {
        axisLine: {
          onZero: false,
        },
        position: 'bottom',
        // min: Math.min.apply(null, array_col) - 0.5,
        // max: Math.max.apply(null, array_col) + 0.5,
      },
      yAxis: {
        axisLine: {
          onZero: false,
        },
        position: 'left',
        min: Math.min.apply(null, array_row) - 5,
        max: Math.max.apply(null, array_row) + 5,
      },
      grid: {
        top: '7%',
        left: '10%',
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
      legend: {
        show: true,
        bottom: '0%',
        data: legendData,
        selectedMode: 'single',
        selected: {
          V1_Adult_Mouse_Brain: true,
          BayesSpace: false,
        },
        textStyle: {
          fontSize: 16,
        },
      },
      series: mySeries,
      // graphic: [
      //   {
      //     type: 'group',
      //     left: '10%',
      //     top: 'center',
      //     children: [
      //       {
      //         type: 'rect',
      //         z: 0,
      //         left: 'center',
      //         top: 'middle',
      //         shape: {
      //           width: 240,
      //           height: 90,
      //         },
      //         style: {
      //           color: 'rgba(0,0,0,0)',
      //           fill: '#fff',
      //           stroke: '#555',
      //           lineWidth: 1,
      //           shadowBlur: 8,
      //           shadowOffsetX: 3,
      //           shadowOffsetY: 3,
      //           shadowColor: 'rgba(0,0,0,0.2)',
      //         },
      //       },
      //       {
      //         type: 'image',
      //         z: 0,
      //         left: 'center',
      //         top: 'middle',
      //         style: {
      //           image:
      //             'https://pic1.zhimg.com/80/v2-884c88e9ddd12fd94418c0359a43c148_1440w.webp',

      //           width: 220,
      //         },
      //       },
      //     ],
      //   },
      // ],
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

SpScatter.defaultProps = {
  theme: 'dark',
  title: 'SpScatter',
  height: '35rem',
  width: '35rem',
  margin: '2rem',
}

SpScatter.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default SpScatter
