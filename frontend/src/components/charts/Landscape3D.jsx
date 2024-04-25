import React, { useRef, useEffect, useState } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import { TooltipComponent, VisualMapComponent } from 'echarts/components'
import { CanvasRenderer } from 'echarts/renderers'
import { Scatter3DChart } from 'echarts-gl/charts'
import { Grid3DComponent } from 'echarts-gl/components'
import { Button, Stack, Badge, Form, Accordion } from 'react-bootstrap'

echarts.use([
  TooltipComponent,
  VisualMapComponent,
  Grid3DComponent,
  Scatter3DChart,
  CanvasRenderer,
])

export const Landscape3D = ({
  theme,
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
    var myChart = echarts.init(chartRef.current, theme) //init the echart container
    var option
    var data
    var config = (app.config = {
      xAxis3D: 'X',
      yAxis3D: 'Y',
      zAxis3D: 'Z',
      color: 'color',
      cell_type: 'cell_type',
      //symbolSize: 'vitaminc',
      onChange: function () {
        if (data) {
          myChart.setOption({
            xAxis3D: {
              name: config.xAxis3D,
            },
            yAxis3D: {
              name: config.yAxis3D,
            },
            zAxis3D: {
              name: config.zAxis3D,
            },
            series: {
              dimensions: [
                config.xAxis3D,
                config.yAxis3D,
                config.zAxis3D,
                config.color,
                config.cell_type,
              ],
              data: _data.map((item) => {
                return [
                  item[config.xAxis3D],
                  item[config.yAxis3D],
                  item[config.zAxis3D],
                  item[config.color],
                  item[config.cell_type],
                ]
              }),
            },
          })
        }
      },
    })

    var _data = require('../../assets/data/PDAC-3D.json')

    var ct_datas = []

    var currentGroup
    for (var i = 0; i < _data.length; i++) {
      // seperate the _data into subgroups by cell-types
      var item = _data[i]
      if (!currentGroup || currentGroup.cell_type !== item.cell_type) {
        currentGroup = { cell_type: item.cell_type, points: [] }
        ct_datas.push(currentGroup)
      }
      currentGroup.points.push(item)
    }
    var sereisArray = []
    console.log(visible)
    for (var i = 0; i < ct_datas.length; i++) {
      var ct_data = ct_datas[i]
      sereisArray.push({
        show: visible,
        tooltip: {},
        name: ct_data.cell_type,
        type: 'scatter3D',
        dimensions: [
          config.xAxis3D,
          config.yAxis3D,
          config.zAxis3D,
          config.color,
          config.cell_type,
        ],
        data: ct_data.points.map((item) => {
          return [
            item[config.xAxis3D],
            item[config.yAxis3D],
            item[config.zAxis3D],
            item[config.color],
            item[config.cell_type],
          ]
        }),
        symbolSize: visible ? 10 : 30,
        // symbol: 'triangle',a
        itemStyle: {
          color: (params) => {
            return params.data[3]
          },
        },
        emphasis: {
          focus: 'series',
        },
      })
    }
    console.log(sereisArray)
    myChart.setOption({
      tooltip: {},
      xAxis3D: {
        name: config.xAxis3D,
        type: 'value',
      },
      yAxis3D: {
        name: config.yAxis3D,
        type: 'value',
      },
      zAxis3D: {
        name: config.zAxis3D,
        type: 'value',
      },
      grid3D: {
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
        viewControl: {
          autoRotate: true,
          // projection: 'orthographic'
        },
      },
      legend: {
        // add legend components
        show: true,
        data: [
          'Acinar cells',
          'Cancer clone A',
          'Cancer clone B',
          'Fibroblasts',
        ],
        textStyle: {
          color: 'white',
          fontSize: '24',
        },
        orient: 'horizontal',
        left: 'center',
        top: 'bottom',
      },

      toolbox: {
        show: true,
        orient: 'vertical',
        left: 'right',
        top: 'center',
        feature: {
          mark: { show: true },
          dataView: { show: true, readOnly: false },
          magicType: { show: true, type: ['line', 'bar', 'stack'] },
          restore: { show: true },
          saveAsImage: { show: true },
        },
      },
      series: sereisArray,
    })

    option && myChart.setOption(option)

    // 定义一个函数来隐藏和显示系列

    // function toggleSeries(seriesName, visible) {
    //   var action = {
    //     type: visible ? 'hide' : 'show', // 使用 'show' 来显示系列，使用 'hide' 来隐藏系列
    //     seriesIndex: function (idx) {
    //       return (
    //         idx ===
    //         myChart.getOption().series.findIndex(function (series) {
    //           return series.name === seriesName
    //         })
    //       )
    //     },
    //   }

    //   myChart.dispatchAction(action)
    // }
  }, [theme, title, visible])

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

Landscape3D.defaultProps = {
  theme: 'dark',
  title: 'Landscape3D',
  height: '35rem',
  width: '35rem',
  margin: '2rem',
}

Landscape3D.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default Landscape3D
