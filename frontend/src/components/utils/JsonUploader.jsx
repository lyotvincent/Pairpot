import React, { useState, useEffect } from 'react'
import { Upload, Button } from 'antd'
import { UploadOutlined } from '@ant-design/icons'
import axios from 'axios'

const JsonUploader = () => {
  const beforeUpload = (file) => {
    const isJson =
      file.type === 'application/json' || file.name.endsWith('.json')
    if (!isJson) {
      alert('Please select a .json file~')
      return false
    } else {
      return true
    }
  }

  const onUpload = (info) => {
    const file = info.file
    const reader = new FileReader()
    reader.onload = (event) => {
      const jsonData = JSON.parse(event.target.result)
      axios
        .post('http://localhost:5001/upload', {
          data: jsonData,
          headers: {
            'Content-Type': 'application/json',
          },
        })
        .then((response) => console.log(response))
        .catch((error) => console.log(error))
    }
    reader.readAsText(file)
  }

  const onChange = (info) => {
    console.log(info)
  }
  const props = {
    maxCount: 1,
    beforeUpload: beforeUpload,
    customRequest: onUpload,
    onChange: onChange,
    listType: 'text',
  }
  return (
    <Upload {...props}>
      <Button type="primary" block icon={<UploadOutlined />}>
        Upload
      </Button>
    </Upload>
  )
}

export default JsonUploader
