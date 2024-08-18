// our SearchEngine
import React, { useState } from 'react'
import { UserOutlined } from '@ant-design/icons'
import { AutoComplete, Input, Col, Row, Popover, ConfigProvider } from 'antd'
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
const options = [
  {
    label: renderTitle('Spatial Transcriptomics Studies'),
    options: [renderItem('Spatial transcriptomics map of the embryonic mouse brain – a tool to explore neurogenesis', 'STDS0000235'),
    renderItem('Spatial transcriptomics combined with single-cell RNA-sequencing unravels the complex inflammatory cell network in atopic dermatitis [ST]', 'STDS0000212'),
    // renderItem("Spatial transcriptomics of de novo NEPC and ARPC", 'STDS0000227'),
    // renderItem("Spatial resolution of cellular senescence dynamics in colorectal liver metastasis", 'STDS0000219'),
    // renderItem("Spatial transcriptomics of adenoid cystic carcinoma of the lacrimal gland", 'STDS0000221'),
    // renderItem("Spatial transcriptomics (Visium, 10x Genomics) data of Duchenne mouse models", 'STDS0000204'),
    // renderItem("Discovering Haematoma-Stimulated Circuits for Secondary Brain Injury after Intraventricular Haemorrhage by Spatial Transcriptome Analysis", 'STDS0000201'),
    // renderItem("Moxibustion Improves Hypothalamus Aqp4 Polarization in APP/PS1 Mice: Evidence from Spatial Transcriptomics", 'STDS0000189')
    ],
  },
]


const Search = ({ onSearchComplete }) => {
  const [loadText, setLoadText] = useState("Search")
  const navigate = useNavigate()
  const [loading, setLoading] = useState([])
  const [open, setOpen] = useState(false)
  const [status, setStatus] = useState("Success")
  const fetch = useMutation({
    mutationKey:['search_key'],
    mutationFn:(e) =>{
      axios({
        method: 'GET',
        url: '/api/search_key',
        params: {
          type: e.type,
          content: e.content
        },
      }).then((response) => {
        let datas = response.data.data
        onSearchComplete(datas)  // 将数据传递给父组件
        quitLoading(0,setLoading)
        setLoadText("Search")  // 重置按钮文本
        return datas
      }).catch(() => {
        setLoading(0,setLoading)  // 搜索失败时也恢复状态
        setLoadText("Search")
      })
    }
  })
  return (
    <ConfigProvider theme={{
      components: {
        Popover: {
          colorBgElevated: '#ddf7d4',
        }
      }
    }}>
      <AutoComplete
        popupClassName="certain-category-search-dropdown"
        popupMatchSelectWidth="80%"
        style={{
          width: "100%"
        }}
        options={options}
        open={open}
        onBlur={() => { setOpen(false) }}
        onSelect={() => { setOpen(false) }}
        size="large"
      >
        <Input.Search
          size="large"
          placeholder="Search for a Study to Browse."
          enterButton={loadText}
          status={status}
          loading={loading[0]}
          onClick={() => { setOpen(!open) }}
          onSearch={(e) => {
            if(e.length > 0){
                // 必要的状态设置
                setOpen(false)
                setStatus("Success")
                enterLoading(0, setLoading)
                setLoadText("Searching...")

                // 需要在这判断一下搜索的类型
                var res;
                // case1: id搜索
                if(e.startsWith("STDS") || e.startsWith("SCDS")){
                    res={
                        type: 'id',
                        content: e
                    }
                }
                else{ // case2: keyword搜索
                    res={
                        type: 'key',
                        content: e
                    }
                }


                fetch.mutate(res)
            } 
            else{
                setStatus("error")
            }
            
          }} />
      </AutoComplete>
    </ConfigProvider>
  )
}

export default Search
