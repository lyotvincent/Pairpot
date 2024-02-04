import React, { useRef, useEffect } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'

const HistPlot = ({ theme, title, height, width, margin }) => {
  const chartRef = useRef() // get current DOM container
  useEffect(() => {
    var chart = echarts.init(chartRef.current, theme) //init the echart container
    var app = {}
    var option
    const posList = [
      'left',
      'right',
      'top',
      'bottom',
      'inside',
      'insideTop',
      'insideLeft',
      'insideRight',
      'insideBottom',
      'insideTopLeft',
      'insideTopRight',
      'insideBottomLeft',
      'insideBottomRight',
    ]
    app.configParameters = {
      rotate: {
        min: -90,
        max: 90,
      },
      align: {
        options: {
          left: 'left',
          center: 'center',
          right: 'right',
        },
      },
      verticalAlign: {
        options: {
          top: 'top',
          middle: 'middle',
          bottom: 'bottom',
        },
      },
      position: {
        options: posList.reduce(function (map, pos) {
          map[pos] = pos
          return map
        }, {}),
      },
      distance: {
        min: 0,
        max: 100,
      },
    }
    app.config = {
      rotate: 90,
      align: 'left',
      verticalAlign: 'middle',
      position: 'insideBottom',
      distance: 15,
      onChange: function () {
        const labelOption = {
          rotate: app.config.rotate,
          align: app.config.align,
          verticalAlign: app.config.verticalAlign,
          position: app.config.position,
          distance: app.config.distance,
        }
        chart.setOption({
          series: [
            {
              label: labelOption,
            },
            {
              label: labelOption,
            },
            {
              label: labelOption,
            },
            {
              label: labelOption,
            },
          ],
        })
      },
    }
    const labelOption = {
      show: true,
      position: app.config.position,
      distance: app.config.distance,
      align: app.config.align,
      verticalAlign: app.config.verticalAlign,
      rotate: app.config.rotate,
      formatter: '{c}  {name|{a}}',
      fontSize: 16,
      rich: {
        name: {},
      },
    }
    option = {
      title: {
        text: title,
        left: 'center',
        top: '5%',
      },
      grid: {
        height: '70%',
        top: '15%',
      },
      tooltip: {
        trigger: 'axis',
        axisPointer: {
          type: 'shadow',
        },
      },
      legend: {
        show: true,
        data: ['Forest', 'Steppe', 'Desert', 'Wetlande'],
        //bottom: '10%',
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
      xAxis: [
        {
          type: 'category',
          axisTick: { show: false },
          data: ['2012', '2013', '2014', '2015', '2016'],
        },
      ],
      yAxis: [
        {
          type: 'value',
        },
      ],
      series: [
        {
          name: 'Forest',
          type: 'bar',
          barGap: 0,
          label: labelOption,
          emphasis: {
            focus: 'series',
          },
          data: [320, 332, 301, 334, 390],
        },
        {
          name: 'Steppe',
          type: 'bar',
          label: labelOption,
          emphasis: {
            focus: 'series',
          },
          data: [220, 182, 191, 234, 290],
        },
        {
          name: 'Desert',
          type: 'bar',
          label: labelOption,
          emphasis: {
            focus: 'series',
          },
          data: [150, 232, 201, 154, 190],
        },
        {
          name: 'Wetlande',
          type: 'bar',
          label: labelOption,
          emphasis: {
            focus: 'series',
          },
          data: [98, 77, 101, 99, 40],
        },
      ],
    }

    option && chart.setOption(option)
  }, [theme, title])

  return (
    <div
      ref={chartRef}
      className="chart"
      //the target DOM container needs height and width
      style={{ height: height, width: width, margin: margin }}></div>
  )
}

HistPlot.defaultProps = {
  theme: 'dark',
  title: 'HistPlot',
  height: '35rem',
  width: '35rem',
  margin: '2rem',
}

HistPlot.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default HistPlot
