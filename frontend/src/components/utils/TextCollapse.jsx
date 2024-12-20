import React, { useState } from 'react'
import { Button } from 'antd'

const TextCollapse = ({ text, header, threshold = 100 }) => {
  const [detail, setDetail] = useState(false)

  return (
    <div>
      {header}
      {detail ? text : text?.substring(0, threshold)}
      <Button
        onClick={() => setDetail(!detail)}
        type="link"
        size="small"
        theme={{
          components: {
            Button: {
              /* here is your component tokens */
              paddingInline: 0,
            },
          },
        }}>
        {detail ? 'show less' : '...show details'}
      </Button>
    </div>
  )
}

export default TextCollapse
