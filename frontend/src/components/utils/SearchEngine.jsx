// our SearchEngine
import React, { useState } from 'react'
import { UserOutlined } from '@ant-design/icons'
import { AutoComplete, Input, Col, Row, Popover, ConfigProvider, Button } from 'antd'
import { useNavigate, Link } from 'react-router-dom'
import Loading from '../charts/Loading'
import axios from 'axios'
import TextCollapse from './TextCollapse'
import { useMutation } from 'react-query'
const { enterLoading, quitLoading } = Loading
const renderTitle = (title) => (
  <span>
    {title}
    <Link
      style={{
        float: 'right',
      }}
      to="/database"
    >
      more
    </Link>
  </span>
)
const renderItem = (title, key) => ({
  value: key,
  label: (
    <div
      style={{
        display: 'flex',
        justifyContent: 'space-between',
      }}
    >
      <Popover content={title} zIndex={2000}>
        {title.length > 90 ? `${title.substring(0, 90)}...` : title}
      </Popover>
      <span>
        <UserOutlined /> {key}
      </span>
    </div>
  ),
})
const initialOptions = [
  {
    label: renderTitle('Spatial Transcriptomics Studies'),
    options: [renderItem('Spatial transcriptomics map of the embryonic mouse brain – a tool to explore neurogenesis', 'STDS0000235'),
    renderItem('Spatial transcriptomics combined with single-cell RNA-sequencing unravels the complex inflammatory cell network in atopic dermatitis [ST]', 'STDS0000212'),
    renderItem("Spatial transcriptomics of de novo NEPC and ARPC", 'STDS0000227'),
    renderItem("Spatial resolution of cellular senescence dynamics in colorectal liver metastasis", 'STDS0000219'),
    renderItem("Spatial transcriptomics of adenoid cystic carcinoma of the lacrimal gland", 'STDS0000221'),
    renderItem("Spatial transcriptomics (Visium, 10x Genomics) data of Duchenne mouse models", 'STDS0000204'),
    renderItem("Discovering Haematoma-Stimulated Circuits for Secondary Brain Injury after Intraventricular Haemorrhage by Spatial Transcriptome Analysis", 'STDS0000201'),
    renderItem("Moxibustion Improves Hypothalamus Aqp4 Polarization in APP/PS1 Mice: Evidence from Spatial Transcriptomics", 'STDS0000189')
    ],
  },
]


const Search = ({ onSearchComplete, onChange, value }) => { // 增加输入框的value和改变函数onchange
  // const [value, setValue] = useState('')  // 这里定义 value 和 setValue
  const [loadText, setLoadText] = useState("Search")
  const navigate = useNavigate()
  const [loading, setLoading] = useState([])
  const [open, setOpen] = useState(false)
  const [status, setStatus] = useState("Success")
  const [options, setOptions] = useState(initialOptions)

  const fetch = useMutation({
    mutationKey: ['search_key'],
    mutationFn: (e) => {
      console.log(e)
      axios({
        method: 'GET',
        url: '/api/search_key',
        params: {
          type: e.type,
          content: e.content
        },
      }).then((response) => {
        // let datas = response.data.data
        let datas = response.data
        console.log(datas)
        onSearchComplete(datas, e)  // to parent
        quitLoading(0, setLoading)
        setLoadText("Search")
        return datas
      }).catch(() => {
        setLoading(0, setLoading)
        setLoadText("Search")
      })
    }
  })

  // 下拉框根据输入来动态更新
  const handleSearch = (value) => {
    onChange(value) // 更新搜索框状态
    if (value.length > 0) {
        const filteredOptions = initialOptions.map(optionGroup => ({
          ...optionGroup,
          options: optionGroup.options.filter(item => item.label.props.children[0].props.children.toLowerCase().includes(value.toLowerCase()))
        }))
        setOptions(filteredOptions)
    } else {
      setOptions(initialOptions)  // 如果没有输入内容，则重置为初始的选项
    }
  }

  // 检查字符串是不是纯数字
  const isNumeric = (str) => {
    return str.split('').every(char => char >= '0' && char <= '9');
  }


  return (
    <ConfigProvider theme={{
      components: {
        Popover: {
          colorBgElevated: '#ddf7d4',
        }
      }
    }}>
      <Row gutter={8}>
        <Col flex="auto">
          <AutoComplete
            popupClassName="certain-category-search-dropdown"
            popupMatchSelectWidth="80%"
            style={{
              width: "100%"
            }}
            options={options}
            open={open}
            onBlur={() => { setOpen(false) }}
            // onSelect={() => { setOpen(false) }}
            onSelect={(selectedValue) => {  // 点击选项时触发
              onChange(selectedValue)  // 更新父组件的value
              setOpen(false)  // 关闭下拉框
            }}
            onSearch={handleSearch}  // 输入变化
            size="large"
            value={value} // 父组件Database传来的value值
            // onChange={(newValue) => setValue(newValue)} // 更新输入框中的value
          >

            <Input.Search
              size="large"
              placeholder="e.g., STDS0000235, mouse brain, cancer, etc."
              enterButton={loadText}
              status={status}
              loading={loading[0]}
              onClick={() => { setOpen(!open) }}
              onSearch={(e) => {
                if (e.length > 0) {
                  // status settings
                  setOpen(false)
                  setStatus("Success")
                  enterLoading(0, setLoading)
                  setLoadText("Searching...")

                  // query type
                  var res
                  // case1: id
                  if (e.startsWith("STDS") || e.startsWith("SCDS")) {
                    res = {
                      type: 'id',
                      content: e
                    }
                  }
                  else if(isNumeric(e)){
                    // case2: num
                    // 包装成id
                    res = {
                      type: 'num',
                      content: e
                    }
                  }
                  else { 
                    // case3: keyword
                    res = {
                      type: 'key',
                      content: e
                    }
                  }



                  fetch.mutate(res)
                }
                else {
                  setStatus("error")
                }

              }} />

          </AutoComplete>
        </Col>

        <Col>
          <Button
            size="large"
            type="primary"
            loading={loading[0]}
            onClick={() => {
              // onChange('') // 清空搜索框 不需要，父组件已经处理过了！！！
              // setValue('')
              var res = {
                type: 'all',
                content: 'all'
              }
              fetch.mutate(res)
            }}
          >
            Display All
          </Button>
        </Col>
      </Row>

    </ConfigProvider>
  )
}

export default Search
