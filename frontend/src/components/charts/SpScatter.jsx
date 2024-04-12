import React, { useRef, useEffect, useState } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import { TooltipComponent, VisualMapComponent } from 'echarts/components'
import {
  ConfigProvider,
  Button,
  Card,
  Space,
  Popconfirm,
  Input,
  Select,
  notification,
  Popover,
  Upload,
  Segmented,
  InputNumber,
  Slider,
  Switch,
  Col,
  Row,
  Form,
  Flex
} from 'antd'
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
  const [api, contextHolder] = notification.useNotification()
  const [isInit, setInit] = useState(false) // whether echart object is inited
  const [xInv, setxInv] = useState(false)
  const [yInv, setyInv] = useState(false)
  const itemGroupRef = useRef([])

  const setAxis = (
    source,
    dims,
    xName0,
    yName0,
    xName1 = null,
    yName1 = null
  ) => {
    // the axis of left chart
    let xIdx0 = dims.indexOf(xName0)
    let yIdx0 = dims.indexOf(yName0)
    let embd0_x = source.map((item) => {
      return item[xIdx0]
    })
    let embd0_y = source.map((item) => {
      return item[yIdx0]
    })

    // the axis of right chart. By default it's as same as the axis in left chart.
    let xIdx1 = xIdx0
    let yIdx1 = yIdx0
    let embd1_x = embd0_x
    let embd1_y = embd0_y
    if (xName1 !== null && yName1 !== null) {
      xIdx1 = dims.indexOf(xName1)
      yIdx1 = dims.indexOf(yName1)
      embd1_x = source.map((item) => {
        return item[xIdx1]
      })
      embd1_y = source.map((item) => {
        return item[yIdx1]
      })
    }

    return {
      xAxis: [
        {
          id: 0,
          gridIndex: 0,
          name: xName0,
          nameLocation: 'middle',
          nameGap: 23,
          nameTextStyle: {
            fontSize: 16,
          },
          axisLine: {
            onZero: false,
          },
          axisLabel: {
            fontSize: 14,
          },
          position: 'bottom',
          inverse: xInv,
          min: Math.ceil(Math.min.apply(null, embd0_x) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd0_x) + 1.5),
        },
        {
          id: 1,
          gridIndex: 1,
          name: xName1,
          nameLocation: 'middle',
          nameGap: 23,
          nameTextStyle: {
            fontSize: 16,
          },
          axisLine: {
            onZero: false,
          },
          axisLabel: {
            fontSize: 14,
          },
          position: 'bottom',
          inverse: xInv,
          min: Math.ceil(Math.min.apply(null, embd1_x) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd1_x) + 1.5),
        },
      ],
      yAxis: [
        {
          id: 0,
          gridIndex: 0,
          name: yName0,
          nameLocation: yInv ? 'start' : 'end',
          nameTextStyle: {
            fontSize: 16,
          },
          axisLine: {
            onZero: false,
          },
          axisLabel: {
            fontSize: 14,
          },
          position: 'left',
          inverse: yInv,
          min: Math.ceil(Math.min.apply(null, embd0_y) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd0_y) + 1.5),
        },
        {
          id: 1,
          gridIndex: 1,
          name: yName1,
          nameLocation: yInv ? 'start' : 'end',
          nameTextStyle: {
            fontSize: 16,
          },
          axisLine: {
            onZero: false,
          },
          axisLabel: {
            fontSize: 14,
          },
          position: 'left',
          inverse: yInv,
          min: Math.ceil(Math.min.apply(null, embd1_y) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd1_y) + 1.5),
        },
      ],
    }
  }

  const setItemGroup = (source, annoIdx, type = 'categories') => {
    let itemGroup = []
    let annotations = []
    if (type === 'categories') {
      annotations = [
        ...new Set(
          source.map((item) => {
            return item[annoIdx]
          })
        ),
      ]

      annotations = annotations.sort()
      for (let anno of annotations) {
        itemGroup.push(
          source
            .filter((item) => {
              return item[annoIdx] === anno
            })
            .map((item) => item[item.length - 1])
        )
      }
    } else {
      itemGroup.push(source.map((item) => item[item.length - 1]))
    }

    itemGroupRef.current = itemGroup
    return annotations
  }

  const setDataset = (source, dims, annoName, annotations) => {
    let _datasets = []
    _datasets.push({
      dimensions: dims,
      source: source,
    })
    if (annotations.length === 0) {
      // if no annotations, all the data are annotated by a label.
      _datasets.push({
        // 这个 dataset 的 index 是 `1`。
        transform: {
          type: 'sort',
          config: { dimension: annoName, order: 'desc' },
        },
      })
    } else {
      for (let anno of annotations) {
        _datasets.push({
          // 这个 dataset 的 index 是 `1`。
          transform: {
            type: 'filter',
            config: { dimension: annoName, value: anno },
          },
        })
      }
    }
    return _datasets
  }

  useEffect(() => {
    if (isInit) {
      console.log("Is Inited.")
      var myChart = echarts.getInstanceByDom(chartRef.current)
    } else {
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
