package kip.nativelib;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;
import com.sun.jna.Native;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import java.util.Collections;
import java.util.Map;

public class FreeableMemory extends Memory {
    
    private static Map<ByteBuffer, FreeableMemory> buffers = Collections.synchronizedMap(new WeakIdentityHashMap<ByteBuffer, FreeableMemory>());
    
    /** Force cleanup of memory that has associated NIO Buffers which have
        been GC'd.
    */ 
    public static void purge() {
        buffers.size();
    }

    /** Free the native memory and set peer to zero */
    public void dispose() {
        super.dispose();
    }

    /**
     * Get a ByteBuffer mapped to a portion of this memory.  
     * We keep a weak reference to all ByteBuffers provided so that this
     * memory object is not GC'd while there are still implicit outstanding
     * references to it (it'd be nice if we could attach our own reference to
     * the ByteBuffer, but the VM generates the object so we have no control
     * over it).
     *
     * @param offset byte offset from pointer to start the buffer
     * @param length Length of ByteBuffer
     * @return a direct ByteBuffer that accesses the memory being pointed to, 
     */
    public ByteBuffer getByteBuffer(long offset, long length) {
        boundsCheck(offset, length);
        ByteBuffer b = Native.getDirectByteBuffer(peer + offset, length).order(ByteOrder.nativeOrder());
        // Ensure this Memory object will not be GC'd (and its memory freed)
        // if the Buffer is still extant.
        buffers.put(b, this);
        return b;
    }
}