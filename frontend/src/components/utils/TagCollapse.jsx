import React, { useState } from 'react'
import { Button } from 'antd'

const TagCollapse = ({ tags, threshold = 5, prefix = '' }) => {
  const [detail, setDetail] = useState(false)

  return (
    <div>
      {prefix}
      {detail ? tags : tags.slice(0, threshold)}
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
        {detail ? '<< show less' : 'show more >>'}
      </Button>
    </div>
  )
}

export default TagCollapse
