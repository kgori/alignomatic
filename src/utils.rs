use anyhow::{anyhow, Result};
use std::fs::File;
use std::os::unix::io::AsRawFd;
pub fn silence_stderr<T, F>(f: F) -> Result<T>
where
    F: FnOnce() -> T,
{
    let dev_null = File::open("/dev/null")?;
    let dev_null_fd = dev_null.as_raw_fd();
    let stderr_fd = unsafe { libc::dup(libc::STDERR_FILENO) };

    if stderr_fd < 0 {
        return Err(anyhow!("Failed to duplicate stderr file descriptor"));
    }

    if unsafe { libc::dup2(dev_null_fd, libc::STDERR_FILENO) } < 0 {
        return Err(anyhow!("Failed to redirect stderr to /dev/null"));
    }

    let result = f();

    if unsafe { libc::dup2(stderr_fd, libc::STDERR_FILENO) } < 0 {
        return Err(anyhow!("Failed to restore stderr"));
    }

    Ok(result)
}
